// native/file_io.cpp — reading and writing whole files, with gzip (see file_io.hpp).

#include "file_io.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include <zlib.h>

namespace lgca {

bool read_file(const std::string& path, std::size_t limit, std::vector<unsigned char>& bytes,
               std::string& error) {
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> file(std::fopen(path.c_str(), "rb"),
                                                               std::fclose);
    if (!file) {
        error = std::strerror(errno);
        return false;
    }
    bytes.clear();
    std::vector<unsigned char> in(1 << 20);
    std::size_t have = std::fread(in.data(), 1, in.size(), file.get());
    const bool gzip = have >= 3 && in[0] == 0x1f && in[1] == 0x8b && in[2] == 8;
    const std::string too_large =
        "more than " + std::to_string(limit) + " bytes" + (gzip ? " after decompression" : "");
    const std::size_t gzip_limit = limit + limit / 1024 + (1 << 20);
    std::size_t gzip_read = have;
    // The next chunk of the file into `in`; false on a read error or a gzip file too long.
    const auto read_more = [&]() {
        have = std::fread(in.data(), 1, in.size(), file.get());
        if (std::ferror(file.get())) {
            error = "read error";
            return false;
        }
        gzip_read += have;
        if (gzip && gzip_read > gzip_limit) {
            error = "the gzip file is longer than " + std::to_string(gzip_limit) + " bytes";
            return false;
        }
        return true;
    };
    if (std::ferror(file.get())) {
        error = "read error";
        return false;
    }

    if (!gzip) {
        while (have > 0) {
            if (have > limit - bytes.size()) {
                error = too_large;
                return false;
            }
            bytes.insert(bytes.end(), in.begin(), in.begin() + have);
            if (!read_more()) return false;
        }
        return true;
    }

    z_stream z{};
    if (inflateInit2(&z, 16 + MAX_WBITS) != Z_OK) {   // gzip format only
        error = "zlib could not start";
        return false;
    }
    const std::unique_ptr<z_stream, int (*)(z_stream*)> end(&z, inflateEnd);
    std::vector<unsigned char> out(1 << 20);
    z.next_in = in.data();
    z.avail_in = static_cast<uInt>(have);
    while (true) {
        if (z.avail_in == 0) {
            if (!read_more()) return false;
            if (have == 0) {
                error = "the gzip stream is cut off (unexpected end of file)";
                return false;
            }
            z.next_in = in.data();
            z.avail_in = static_cast<uInt>(have);
        }
        z.next_out = out.data();
        z.avail_out = static_cast<uInt>(out.size());
        const int code = inflate(&z, Z_NO_FLUSH);
        if (code != Z_OK && code != Z_STREAM_END && code != Z_BUF_ERROR) {
            error = std::string("damaged gzip stream (") + (z.msg ? z.msg : "zlib error") + ")";
            return false;
        }
        const std::size_t produced = out.size() - z.avail_out;
        if (produced > limit - bytes.size()) {
            error = too_large;
            return false;
        }
        bytes.insert(bytes.end(), out.begin(), out.begin() + produced);
        if (code != Z_STREAM_END) continue;
        // The end of one gzip member. What follows decides, as in Node: nothing (done),
        // a zero byte (padding: the rest is ignored), or anything else, read as the next
        // member, so that data which is not a gzip member is an error.
        if (z.avail_in == 0) {
            if (!read_more()) return false;
            if (have == 0) return true;
            z.next_in = in.data();
            z.avail_in = static_cast<uInt>(have);
        }
        if (*z.next_in == 0) return true;
        inflateReset(&z);
    }
}

bool gzip_compress(const void* data, std::size_t size, std::vector<unsigned char>& out,
                   std::string& error) {
    z_stream z{};
    // windowBits 16 + 15: the gzip format; memLevel 8 and the default strategy are
    // zlib's (and Node's) defaults.
    if (deflateInit2(&z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK) {
        error = "zlib could not start";
        return false;
    }
    const std::unique_ptr<z_stream, int (*)(z_stream*)> end(&z, deflateEnd);
    out.clear();
    const unsigned char* in = static_cast<const unsigned char*>(data);
    std::size_t left = size;
    std::vector<unsigned char> chunk(1 << 20);
    int code;
    do {
        // zlib counts in 32-bit units: feed the input in pieces of at most 1 GiB.
        if (z.avail_in == 0 && left > 0) {
            const std::size_t piece = left < (std::size_t(1) << 30) ? left : (std::size_t(1) << 30);
            z.next_in = const_cast<unsigned char*>(in);
            z.avail_in = static_cast<uInt>(piece);
            in += piece;
            left -= piece;
        }
        z.next_out = chunk.data();
        z.avail_out = static_cast<uInt>(chunk.size());
        code = deflate(&z, left == 0 ? Z_FINISH : Z_NO_FLUSH);
        if (code == Z_STREAM_ERROR) {
            error = "zlib failed to compress";
            return false;
        }
        out.insert(out.end(), chunk.data(), chunk.data() + (chunk.size() - z.avail_out));
    } while (code != Z_STREAM_END);
    return true;
}

namespace {

// fsync, tried again when a signal interrupts it.
int fsync_retrying(int fd) {
    int code;
    do code = ::fsync(fd);
    while (code != 0 && errno == EINTR);
    return code;
}

}  // namespace

bool write_file_atomic(const std::string& path, const void* data, std::size_t size,
                       std::string& error) {
    // A new temporary file: the process id keeps two programs apart, the counter two
    // writes of this one. A name already taken (left by a stopped program whose process
    // id is ours now) is passed over. Mode 0666 less the umask, as fopen creates files.
    static unsigned long counter = 0;
    std::string temp;
    int fd = -1;
    for (int attempt = 0; fd < 0; ++attempt) {
        temp = path + ".tmp." + std::to_string(::getpid()) + "." + std::to_string(counter++);
        fd = ::open(temp.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0666);
        if (fd < 0 && (errno != EEXIST || attempt == 99)) {
            error = temp + ": " + std::strerror(errno);
            return false;
        }
    }
    std::FILE* f = ::fdopen(fd, "wb");
    if (!f) {
        error = temp + ": " + std::strerror(errno);
        ::close(fd);
        std::remove(temp.c_str());
        return false;
    }
    // The first step that fails gives the error message.
    const char* failed = nullptr;
    int code = 0;
    const auto step = [&](bool ok, const char* what) {
        if (!ok && !failed) {
            failed = what;
            code = errno;
        }
    };
    step(std::fwrite(data, 1, size, f) == size, "cannot write");
    step(!failed && std::fflush(f) == 0, "cannot write");
    step(!failed && fsync_retrying(fileno(f)) == 0, "cannot flush to the disk");
    step(std::fclose(f) == 0, "cannot close");
    if (failed) {
        error = temp + ": " + failed + " (" + std::strerror(code) + ")";
        std::remove(temp.c_str());
        return false;
    }
    if (std::rename(temp.c_str(), path.c_str()) != 0) {
        error = path + ": " + std::strerror(errno);
        std::remove(temp.c_str());
        return false;
    }
    // The rename is an entry in the folder: flush the folder too, so that after a power
    // cut a file written later (run.json) never names this one while the disk lost it.
    // A file system that cannot flush a folder says so with EINVAL or ENOTSUP; then the
    // rename is on the disk when the system gets to it, which is no error. Any other
    // failure is: the file has its new content, but the caller must not go on as if the
    // disk had it.
    const std::size_t slash = path.find_last_of('/');
    const std::string folder = slash == std::string::npos ? "." : path.substr(0, slash ? slash : 1);
    const int dir_fd = ::open(folder.c_str(), O_RDONLY | O_DIRECTORY);
    if (dir_fd < 0) {
        error = folder + ": cannot open the folder to flush it (" + std::strerror(errno) + ")";
        return false;
    }
    const bool flushed = fsync_retrying(dir_fd) == 0;
    const int why = errno;
    ::close(dir_fd);
    if (!flushed && why != EINVAL && why != ENOTSUP && why != EOPNOTSUPP) {
        error = folder + ": cannot flush the folder to the disk (" + std::strerror(why) + ")";
        return false;
    }
    return true;
}

std::string atomic_temp_target(const std::string& name) {
    // A run of one or more digits ending at `end`; returns where it starts, or npos.
    const auto digits_before = [&name](std::size_t end) {
        std::size_t p = end;
        while (p > 0 && name[p - 1] >= '0' && name[p - 1] <= '9') --p;
        return p < end ? p : std::string::npos;
    };
    std::size_t end = name.size();
    const std::size_t counter = digits_before(end);
    if (counter != std::string::npos && counter > 0 && name[counter - 1] == '.') {
        const std::size_t pid = digits_before(counter - 1);
        if (pid != std::string::npos && pid > 0 && name[pid - 1] == '.') end = pid - 1;
    }
    const std::string suffix = ".tmp";
    if (end <= suffix.size() || name.compare(end - suffix.size(), suffix.size(), suffix) != 0) return "";
    return name.substr(0, end - suffix.size());
}

std::vector<std::string> folder_entries(const std::string& dir) {
    namespace fs = std::filesystem;
    std::vector<std::string> names;
    std::error_code ec;
    for (fs::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
        names.push_back(it->path().filename().string());
    if (ec) throw std::runtime_error(dir + ": cannot list the folder: " + ec.message());
    std::sort(names.begin(), names.end());
    return names;
}

FolderLock::FolderLock(const std::string& dir, const std::string& busy) {
    fd_ = ::open(dir.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (fd_ < 0)
        throw std::runtime_error(dir + ": cannot open the folder to lock it (" + std::strerror(errno) +
                                 ")");
    int code;
    do code = ::flock(fd_, LOCK_EX | LOCK_NB);
    while (code != 0 && errno == EINTR);
    if (code != 0) {
        const int why = errno;
        ::close(fd_);
        if (why == EWOULDBLOCK) throw std::runtime_error(busy);
        throw std::runtime_error(dir + ": cannot lock the folder (" + std::strerror(why) + ")");
    }
}

FolderLock::~FolderLock() { ::close(fd_); }

}  // namespace lgca
