/*
 * Minimal implementation of some terminal colors.
 *
 * Inspired by the termcolor library by Ihor Kalnytskyi.
 */

#include <iostream>
#include <string_view>
#include <unistd.h>

namespace termcolor {
    using std::ostream;

    // Send a control sequence only if writing to a tty.
    static ostream& _send(ostream& stream, const char* sequence) {
        bool is_a_tty = (&stream == &std::cout && isatty(STDOUT_FILENO))
                     || (&stream == &std::cerr && isatty(STDERR_FILENO));
        if (is_a_tty)
            stream << std::string_view(sequence);
        return stream;
    }

    static ostream& reset(ostream& stream)      { return _send(stream, "\x1b[0m"); }
    static ostream& blink(ostream& stream)      { return _send(stream, "\x1b[5m"); }
    static ostream& bright_red(ostream& stream) { return _send(stream, "\x1b[91m"); }
};
