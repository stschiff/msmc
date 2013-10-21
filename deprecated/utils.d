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
 
import std.json;
import std.stdio;
import std.math;
import std.traits;
import std.exception;
import std.conv;
import std.typecons;
import std.regex;
import std.algorithm;
import std.array;
import model.msmc_model;
import model.time_intervals;

JSONValue makeJSON(T)(T value) if(isIntegral!T) {
  JSONValue json;
  json.type = JSON_TYPE.INTEGER;
  json.integer = value;
  return json;
}

unittest {
  writeln("testing makeJSON, integral type");
  auto json = makeJSON(-4);
  assert(json.integer == -4);
  json = makeJSON(3U);
  assert(json.integer == 3);
}

JSONValue makeJSON(T)(T value) if(isFloatingPoint!T) {
  enforce(!isNaN(value), "makeJSON!double can't handle NaN");
  JSONValue json;
  json.type = JSON_TYPE.FLOAT;
  if(isInfinity(value)) {
    json.floating = value > 0.0 ? double.max : -double.max;
  }
  else {
    json.floating = value;
  }
  return json;
}

unittest {
  writeln("testing makeJSON(float)");
  auto json = makeJSON(4.4f);
  assert(approxEqual(json.floating, 4.4, 1.0e-6, 0.0));
  json = makeJSON(4.4);
  assert(approxEqual(json.floating, 4.4, 1.0e-8, 0.0));
}

JSONValue makeJSON(T)(T value) if(isSomeString!T) {
  JSONValue json;
  json.type = JSON_TYPE.STRING;
  json.str = value.idup;
  return json;
}

unittest {
  writeln("testing makeJSON(string)");
  auto json = makeJSON("hello");
  assert(json.str == "hello");
}

JSONValue makeJSON(T)(T value) if(is(T == bool)) {
  JSONValue json;
  json.type = value ? JSON_TYPE.TRUE : JSON_TYPE.FALSE;
  return json;
}

unittest {
  writeln("testing makeJSON(bool)");
  auto json = makeJSON(true);
  assert(json.type == JSON_TYPE.TRUE);
}

JSONValue makeJSON(T)(T array) if(isArray!T && !isSomeString!T) {
  JSONValue json;
  json.type = JSON_TYPE.ARRAY;
  json.array = new JSONValue[array.length];
  foreach(i, value; array) {
    json[i] = makeJSON(value);
  }
  return json;
}

unittest {
  writeln("testing makeJSON(double[])");
  auto json = makeJSON([1, 2, 3]);
  assert(json[0].integer == 1);
  assert(json[1].integer == 2);
  json = makeJSON([1.1, 2.2, 3.3]);
  assert(approxEqual(json[0].floating, 1.1, 1.0e-8, 0.0));
  assert(approxEqual(json[1].floating, 2.2, 1.0e-8, 0.0));
  json = makeJSON([[1, 2], [3, 4]]);
  assert(json[0][0].integer == 1);
  assert(json[0][1].integer == 2);
  assert(json[1][0].integer == 3);
  assert(json[1][1].integer == 4);

}

JSONValue makeJSON(T)(T c) if(is(T == char)) {
  return makeJSON([c]);
}

unittest {
  writeln("testing makeJSON(char)");
  auto json = makeJSON('h');
  assert(json.str == "h");
}

T fromJSON(T)(JSONValue json) if(isFloatingPoint!T) {
  T ret;
  enforce(json.type == JSON_TYPE.FLOAT || json.type == JSON_TYPE.INTEGER, text("fromJSON: found type ", json.type, ", but need FLOAT or INTEGER"));
  if(json.type == JSON_TYPE.FLOAT)
    ret = json.floating;
  else
    ret = json.integer;
  return ret;
}

unittest {
  writeln("testing fromJSON!double");
  auto json = makeJSON(4.56);
  assert(fromJSON!double(json) == 4.56);
  assert(fromJSON!float(json) == 4.56f);
  json = makeJSON(5);
  assert(fromJSON!double(json) == 5.0);
  assert(fromJSON!float(json) == 5.0f);
}

T fromJSON(T)(JSONValue json) if(isIntegral!T) {
  enforce(json.type == JSON_TYPE.INTEGER);
  return to!T(json.integer);
}

unittest {
  writeln("testing fromJSON!int");
  auto json = makeJSON(6);
  assert(fromJSON!uint(json) == 6U);
}

T fromJSON(T)(JSONValue json) if(isSomeString!T) {
  enforce(json.type == JSON_TYPE.STRING);
  return to!T(json.str);
}

unittest {
  writeln("testing fromJSON!string");
  auto json = makeJSON("hello");
  assert(fromJSON!string(json) == "hello");
}

T fromJSON(T)(JSONValue json) if(isArray!T && !isSomeString!T) {
  T ret;
  assert(json.type == JSON_TYPE.ARRAY, text(json.type));
  foreach(val; json.array)
    ret ~= fromJSON!(ForeachType!T)(val);
  return ret;
}

unittest {
  writeln("testing fromJSON!(double[])");
  auto json = makeJSON([2.2, 4.4, 5.5]);
  
  auto array = fromJSON!(double[])(json);
  assert(array == [2.2, 4.4, 5.5]);
}

unittest {
  writeln("testing fromJSON!(double[][])");
  auto json = makeJSON([[1.0, 2.0], [3.0, 4.0]]);
  auto matrix = fromJSON!(double[][])(json);
  assert(matrix == [[1.0, 2.0], [3.0, 4.0]]);
}

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
