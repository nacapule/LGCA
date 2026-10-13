// native/json.hpp — a small JSON reader, for run.json files.
//
// The `run` mode (run_cli.cpp) reads an existing run.json to resume a job. It needs a
// few values (the job's settings, the last checkpoint) and the exact text of the
// `series` rows and `snaps` entries it keeps, so every value remembers where its text
// starts and ends in the input.
//
// It reads standard JSON (RFC 8259), as JavaScript's JSON.parse does: objects, arrays,
// strings with every escape (\uXXXX included), numbers, true, false, null; whitespace
// between tokens; nothing but whitespace after the value. Differences that do not
// matter for run.json: nesting deeper than 200 levels is refused, and so is a text of
// more than 1,500,000 values; a number too large for a double reads as infinity, as in
// JavaScript; a \u escape of a lone surrogate half is kept as its 3-byte form (it can
// never equal the plain ASCII values compared here).

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace lgca {

struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;
    bool boolean = false;            // Bool
    double number = 0;               // Number
    std::string text;                // String: the decoded string (UTF-8)
    std::vector<JsonValue> items;    // Array: the elements; Object: the member values
    std::vector<std::string> keys;   // Object: the member names, keys[i] for items[i]
    std::size_t begin = 0, end = 0;  // where the value's text is in the input: [begin, end)

    bool is_number() const { return type == Type::Number; }
    bool is_string() const { return type == Type::String; }
    bool is_array() const { return type == Type::Array; }
    bool is_object() const { return type == Type::Object; }

    // An object's member `key`, or null if there is none (or this is not an object).
    // With a name given twice, the last one counts, as in JSON.parse.
    const JsonValue* get(const std::string& key) const;
};

// Parse the whole of `input`. On an error, returns false and describes it in `error`
// (with the byte offset).
bool parse_json(const std::string& input, JsonValue& value, std::string& error);

}  // namespace lgca
