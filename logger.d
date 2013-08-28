import std.stdio;

File logFile;

static this() {
  logFile = File("/dev/null", "w");
}

void logInfo(string str) {
  logFile.write(str);
  logFile.flush();
  stderr.write(str);
}

