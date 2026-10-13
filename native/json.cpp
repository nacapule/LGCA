// native/json.cpp — a small JSON reader (see json.hpp).

#include "json.hpp"

#include <cerrno>
#include <cstdint>
#include <cstdlib>

namespace lgca {

const JsonValue* JsonValue::get(const std::string& key) const {
    if (type != Type::Object) return nullptr;
    for (std::size_t i = keys.size(); i-- > 0;)
        if (keys[i] == key) return &items[i];
    return nullptr;
}

namespace {

constexpr int MAX_DEPTH = 200;
// At most this many values (each object, array, string, number, true, false or null
// counts as one). A run.json to step 2560000, the last checkpoint, holds about 525,000,
// and a trace.json about 623,000; the cap keeps a damaged file of many short values
// ("[0,0,0,...") from taking gigabytes of memory: each value takes about 100 bytes once
// read, so at most about 150 MB.
constexpr std::size_t MAX_VALUES = 1500000;

class Parser {
public:
    explicit Parser(const std::string& input) : s_(input) {}

    bool parse(JsonValue& value, std::string& error) {
        skip_space();
        if (!parse_value(value, 0)) {
            error = error_;
            return false;
        }
        skip_space();
        if (p_ != s_.size()) {
            error = "unexpected text after the JSON value at byte " + std::to_string(p_);
            return false;
        }
        return true;
    }

private:
    const std::string& s_;
    std::size_t p_ = 0;
    std::size_t values_ = 0;  // values begun so far
    std::string error_;

    bool fail(const std::string& what) {
        error_ = what + " at byte " + std::to_string(p_);
        return false;
    }

    // JSON whitespace: space, tab, line feed, carriage return.
    void skip_space() {
        while (p_ < s_.size() &&
               (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' || s_[p_] == '\r'))
            ++p_;
    }

    bool literal(const char* word) {
        std::size_t q = p_;
        for (const char* c = word; *c; ++c, ++q)
            if (q >= s_.size() || s_[q] != *c) return fail("invalid literal");
        p_ = q;
        return true;
    }

    bool parse_value(JsonValue& v, int depth) {
        if (depth > MAX_DEPTH) return fail("nesting too deep");
        if (++values_ > MAX_VALUES) return fail("more than " + std::to_string(MAX_VALUES) + " values");
        if (p_ >= s_.size()) return fail("unexpected end of the input");
        v = JsonValue();
        v.begin = p_;
        const char c = s_[p_];
        bool ok;
        if (c == '{') {
            v.type = JsonValue::Type::Object;
            ok = parse_object(v, depth);
        } else if (c == '[') {
            v.type = JsonValue::Type::Array;
            ok = parse_array(v, depth);
        } else if (c == '"') {
            v.type = JsonValue::Type::String;
            ok = parse_string(v.text);
        } else if (c == 't') {
            v.type = JsonValue::Type::Bool;
            v.boolean = true;
            ok = literal("true");
        } else if (c == 'f') {
            v.type = JsonValue::Type::Bool;
            ok = literal("false");
        } else if (c == 'n') {
            ok = literal("null");
        } else if (c == '-' || (c >= '0' && c <= '9')) {
            v.type = JsonValue::Type::Number;
            ok = parse_number(v.number);
        } else {
            ok = fail("unexpected character");
        }
        v.end = p_;
        return ok;
    }

    bool parse_object(JsonValue& v, int depth) {
        ++p_;   // '{'
        skip_space();
        if (p_ < s_.size() && s_[p_] == '}') {
            ++p_;
            return true;
        }
        while (true) {
            skip_space();
            if (p_ >= s_.size() || s_[p_] != '"') return fail("expected a member name");
            std::string key;
            if (!parse_string(key)) return false;
            skip_space();
            if (p_ >= s_.size() || s_[p_] != ':') return fail("expected ':'");
            ++p_;
            skip_space();
            v.keys.push_back(std::move(key));
            v.items.emplace_back();
            if (!parse_value(v.items.back(), depth + 1)) return false;
            skip_space();
            if (p_ < s_.size() && s_[p_] == ',') {
                ++p_;
                continue;
            }
            if (p_ < s_.size() && s_[p_] == '}') {
                ++p_;
                return true;
            }
            return fail("expected ',' or '}'");
        }
    }

    bool parse_array(JsonValue& v, int depth) {
        ++p_;   // '['
        skip_space();
        if (p_ < s_.size() && s_[p_] == ']') {
            ++p_;
            return true;
        }
        while (true) {
            skip_space();
            v.items.emplace_back();
            if (!parse_value(v.items.back(), depth + 1)) return false;
            skip_space();
            if (p_ < s_.size() && s_[p_] == ',') {
                ++p_;
                continue;
            }
            if (p_ < s_.size() && s_[p_] == ']') {
                ++p_;
                return true;
            }
            return fail("expected ',' or ']'");
        }
    }

    // Four hex digits of a \u escape.
    bool hex4(uint32_t& code) {
        if (s_.size() - p_ < 4) return fail("cut-off \\u escape");
        code = 0;
        for (int q = 0; q < 4; ++q) {
            const char c = s_[p_++];
            code <<= 4;
            if (c >= '0' && c <= '9') code |= c - '0';
            else if (c >= 'a' && c <= 'f') code |= c - 'a' + 10;
            else if (c >= 'A' && c <= 'F') code |= c - 'A' + 10;
            else return fail("invalid \\u escape");
        }
        return true;
    }

    static void append_utf8(std::string& out, uint32_t code) {
        if (code < 0x80) {
            out += static_cast<char>(code);
        } else if (code < 0x800) {
            out += static_cast<char>(0xc0 | (code >> 6));
            out += static_cast<char>(0x80 | (code & 0x3f));
        } else if (code < 0x10000) {
            out += static_cast<char>(0xe0 | (code >> 12));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3f));
            out += static_cast<char>(0x80 | (code & 0x3f));
        } else {
            out += static_cast<char>(0xf0 | (code >> 18));
            out += static_cast<char>(0x80 | ((code >> 12) & 0x3f));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3f));
            out += static_cast<char>(0x80 | (code & 0x3f));
        }
    }

    bool parse_string(std::string& out) {
        ++p_;   // the opening quote
        out.clear();
        while (true) {
            if (p_ >= s_.size()) return fail("unterminated string");
            const unsigned char c = static_cast<unsigned char>(s_[p_]);
            if (c == '"') {
                ++p_;
                return true;
            }
            if (c < 0x20) return fail("control character in a string");
            if (c != '\\') {
                out += static_cast<char>(c);
                ++p_;
                continue;
            }
            ++p_;   // the backslash
            if (p_ >= s_.size()) return fail("unterminated string");
            const char e = s_[p_++];
            switch (e) {
                case '"': out += '"'; break;
                case '\\': out += '\\'; break;
                case '/': out += '/'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    uint32_t code;
                    if (!hex4(code)) return false;
                    // A surrogate pair written as two escapes is one character.
                    if (code >= 0xd800 && code <= 0xdbff && s_.size() - p_ >= 6 &&
                        s_[p_] == '\\' && s_[p_ + 1] == 'u') {
                        const std::size_t save = p_;
                        p_ += 2;
                        uint32_t low;
                        if (!hex4(low)) return false;
                        if (low >= 0xdc00 && low <= 0xdfff)
                            code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                        else
                            p_ = save;   // not a pair: the second escape is read on its own
                    }
                    append_utf8(out, code);
                    break;
                }
                default:
                    --p_;
                    return fail("invalid escape in a string");
            }
        }
    }

    // A number in JSON's grammar: -?(0|[1-9][0-9]*)(.[0-9]+)?([eE][+-]?[0-9]+)?, read
    // with strtod (correctly rounded, as JSON.parse reads it).
    bool parse_number(double& value) {
        const std::size_t start = p_;
        const auto digits = [this]() {
            const std::size_t from = p_;
            while (p_ < s_.size() && s_[p_] >= '0' && s_[p_] <= '9') ++p_;
            return p_ > from;
        };
        if (s_[p_] == '-') ++p_;
        if (p_ < s_.size() && s_[p_] == '0') {
            ++p_;
        } else if (!digits()) {
            return fail("invalid number");
        }
        if (p_ < s_.size() && s_[p_] == '.') {
            ++p_;
            if (!digits()) return fail("invalid number");
        }
        if (p_ < s_.size() && (s_[p_] == 'e' || s_[p_] == 'E')) {
            ++p_;
            if (p_ < s_.size() && (s_[p_] == '+' || s_[p_] == '-')) ++p_;
            if (!digits()) return fail("invalid number");
        }
        const std::string text = s_.substr(start, p_ - start);
        errno = 0;
        value = std::strtod(text.c_str(), nullptr);   // ERANGE: +-infinity or a tiny value, as in JavaScript
        return true;
    }
};

}  // namespace

bool parse_json(const std::string& input, JsonValue& value, std::string& error) {
    Parser parser(input);
    return parser.parse(value, error);
}

}  // namespace lgca
