// Used only on non-Windows builds (see Makefile): the embedded HTTP server relies on WinSock2.
#include "web_server.hpp"
#include <iostream>

namespace genomix {
void WebServer::start(int, const std::string&, const std::string&) {
    std::cerr << "[ERROR] The embedded Web UI is only available in the Windows (WinSock) build.\n";
}
}
