import std.conv;
import std.typecons;
import std.math;
import std.algorithm;
import std.stdio;
import std.regex;
import std.json;
import std.exception;
import std.array;
import msmc_model;
import utils;
import time_intervals;


size_t[] parseCommaSeparatedArray(string arrayString) {
  enforce(match(arrayString, r"^\d+[,\d+]+"), text("illegal array string: ", arrayString));
  auto splitted = std.string.split(arrayString, ",");
  auto mappedNrs = map!"to!size_t(a)"(splitted);
  return array(mappedNrs);
}
  
unittest {
  writeln("test parseCommaSeparatedArray");
  auto array = parseCommaSeparatedArray("0,0,1,1,2,2");
  assert(array == [0, 0, 1, 1, 2, 2]);
  array = parseCommaSeparatedArray("0,0");
  assert(array == [0, 0]);
  foreach(s; ["23", "2,4.5", "0,1,2"])
    assert(collectExceptionMsg(parseTimeSegmentPattern(s) == text("illegal array string: ", s)));
}

size_t[] parseTimeSegmentPattern(string patternString) {
  enforce(match(patternString, r"^\d+\*\d+[\+\d+\*\d+]*"), text("illegal timeSegmentPattern: ", patternString));
  size_t[] pattern;
  foreach(product; std.string.split(patternString, "+")) {
    auto pair = array(map!"to!size_t(a)"(std.string.split(product, "*")));
    foreach(i; 0 .. pair[0]) {
      pattern ~= pair[1];
    }
  }
  return pattern;
}
  
unittest {
  writeln("test parseTimeSegmentPattern");
  assert(parseTimeSegmentPattern("2*3+1*4") == [3, 3, 4]);
  assert(parseTimeSegmentPattern("1*20") == [20]);
  foreach(pattern; ["aaa1+4+13*2+1+10", "1+2+3", ""])
    assert(collectExceptionMsg(parseTimeSegmentPattern(pattern)) == text("illegal timeSegmentPattern: ", pattern));
}
  

int mOver2(size_t M)
in {
  assert(M >= 2);
}
body {
  return cast(int)(M * (M - 1) / 2);
}

unittest {
  writeln("test mOver2");
  assert(mOver2(2) == 1);
  assert(mOver2(4) == 6);
  assert(mOver2(8) == 28);
}

alias int[4] tree_t;
alias Tuple!(double, int, int) state_triple_t;

state_triple_t readTree(string line) {
  double minTmrca = double.infinity;
  int ind, indd;
  foreach(m; match(line, regex(r"\((?P<ind>[0-9]+):(?P<tmrca>[0-9\.]+),(?P<indd>[0-9]+):[0-9\.]+\)", "g"))) {
    auto c = m.captures;
    auto tmrca = 2.0 * to!double(c["tmrca"]);
    if(tmrca < minTmrca) {
      minTmrca = tmrca;
      ind = to!int(c["ind"]);
      indd = to!int(c["indd"]);
      if(ind > indd) swap(ind, indd);
      assert(ind < indd);
    }
  }
  return state_triple_t(minTmrca, ind, indd);
}
