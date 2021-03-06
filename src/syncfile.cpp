#include "syncfile.hpp"

SyncFile::SyncFile(const std::string& path) : _path(path) { _csv_output.open(_path); }
SyncFile::~SyncFile() { _csv_output.close(); };
bool SyncFile::write(const std::string& data) {
  try {
    std::lock_guard<std::mutex> lock(_writerMutex);
    _csv_output << data << std::endl;
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}

bool SyncFile::write(const csv_data& data) {
  try {
    // Get the file lock for this thread
    std::lock_guard<std::mutex> lock(_writerMutex);
    // put data stream into the file
    _csv_output << data << std::endl;
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}
