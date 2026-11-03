#ifndef WEB_SERVER_HPP
#define WEB_SERVER_HPP

#include <string>

namespace genomix {

class WebServer {
public:
    static void start(int port, const std::string& web_dir, const std::string& data_dir);
};

} // namespace genomix

#endif // WEB_SERVER_HPP
