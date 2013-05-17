import std.exception;
import sub_program;

class SubProgramLauncher {
  
  private static SubProgram[string] subProgramMap;
  private string subProgramName;
  private string[] commandLineArgs;
  
  static void runFromCommandLine(string[] args) {
    auto launcher = new CommandLineLauncher(args);
    launcher.tryToRun();
  }

  static void register(SubProgram fac) {
    subProgramMap[fac.commandName] = fac;
  }

  private this(string[] args) {
    if(args.length < 2)
      exitWithHelpMessage();
    subProgramName = args[1];
    commandLineArgs = args[1..$];
  }

  SubProgram[] sortedSubPrograms() {
    return array(sort!"a.commandName < b.commandName"(subProgramMap.values));
  }
  
  private void tryToRun() {
    if(subProgramExists(subProgramName)) {
      auto program = subProgramMap[subProgramName];
      program.runWithArgs(commandLineArgs);
    }
    else {
      writeln("unknown program ", subProgramName);
      exitWithHelpMessage();
    }
  }
  
  bool subProgramExists(string name) {
    if(name in subProgramMap)
      return true;
    else
      return false;
  }
}
