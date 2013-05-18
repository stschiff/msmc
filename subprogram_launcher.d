/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
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
