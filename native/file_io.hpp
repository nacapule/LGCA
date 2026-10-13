// native/file_io.hpp — reading and writing whole files, with gzip; folder listings and
// folder locks.
//
// Used by the `measure` mode (measure_cli.cpp) to read saved lattices, and by the `run`
// and `trace` modes (run_cli.cpp, trace_cli.cpp) to read and write a job's run.json or
// trace.json and its state-t<T>.bin.gz files.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace lgca {

// The whole content of a file, decompressed if it starts with the gzip bytes 1f 8b 08,
// as Node's gunzipSync decompresses it: several gzip members in a row are read as one;
// after a member, a zero byte ends the reading and the rest of the file is ignored,
// whatever it holds; any other byte must begin another valid member, and a damaged or
// cut-off stream is an error. At most `limit` bytes of content: more is an error. The
// file is read in pieces of 1 MiB, and of a gzip file at most limit + limit/1024 + 1 MiB
// bytes are read (room for stored blocks and a file name in the header), so that input
// without end (empty members, a header that never ends) cannot keep the reading going.
// The reading stops at a zero that ends it: the file beyond the piece holding that zero
// is never read, so it counts toward no limit and its content is never looked at. On an
// error, returns false and says why in `error`.
bool read_file(const std::string& path, std::size_t limit, std::vector<unsigned char>& bytes,
               std::string& error);

// `bytes` gzip-compressed as Node's gzipSync does it by default (zlib, compression level
// 6, no file name, modification time 0). The compressed bytes may differ from Node's
// (Node ships its own zlib build); the content they decompress to is the same.
bool gzip_compress(const void* data, std::size_t size, std::vector<unsigned char>& out,
                   std::string& error);

// Replace the file at `path` with `size` bytes, so that the file is at every moment
// either entirely old or entirely new: the bytes go to a temporary file of this call's
// own, path + ".tmp.<process id>.<n>" (created new, so that two programs writing the
// same file never write into one temporary file), which is flushed to the disk and then
// renamed over `path`; then the folder is flushed, so the rename is on the disk too (on
// a file system that cannot flush a folder, the rename reaches the disk when the system
// writes it; a folder that cannot be opened or flushed for any other reason is an error,
// reported after the rename). An interrupted program leaves at most a stray temporary
// file, which atomic_temp_target recognises, so that its caller can remove it.
bool write_file_atomic(const std::string& path, const void* data, std::size_t size,
                       std::string& error);

// If `name` (a file name, without its folder) is a temporary file that write_file_atomic
// writes, NAME.tmp.<digits>.<digits>, or that earlier builds and the density and alpha
// study's JavaScript runner write, NAME.tmp: NAME, the file it was meant to replace.
// Otherwise "".
std::string atomic_temp_target(const std::string& name);

// The names of the entries of `dir` (without the folder), sorted. Throws
// std::runtime_error if the folder cannot be listed.
std::vector<std::string> folder_entries(const std::string& dir);

// An exclusive lock on a folder: flock on the folder itself (so the folder gets no extra
// file), held while this object lives. The folder is opened by its path, so every name
// that leads to it (a link, another spelling) meets the same lock. Throws
// std::runtime_error with the message `busy` if another process holds the lock, or with
// the reason if the folder cannot be opened or locked. A lock taken through another open
// of the same folder in this process also counts as held: lock each folder once.
class FolderLock {
public:
    FolderLock(const std::string& dir, const std::string& busy);
    ~FolderLock();
    FolderLock(const FolderLock&) = delete;
    FolderLock& operator=(const FolderLock&) = delete;

private:
    int fd_ = -1;
};

}  // namespace lgca
