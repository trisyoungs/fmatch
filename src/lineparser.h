/*
	*** Line Parsing Routines
	*** lineparser.h
	Copyright T. Youngs 2011

	This file is part of FMatch3.

	FMatch3 is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	FMatch3 is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with FMatch3.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FMATCH3_LINEPARSER_H
#define FMATCH3_LINEPARSER_H

#include "dnchar.h"
#include "constants.h"
#include "list.h"
#include <fstream>
#include <iostream>
using namespace std;

#define MAXLINELENGTH 1024

// Line Parser
class LineParser
{
	public:
	// Constructor / Destructor
	LineParser(const char *filename);
	LineParser();
	~LineParser();


	/*
	// Source line/file and read options
	*/
	private:
	// Filename of current file target (if any)
	Dnchar filename_;
	// Line to parse
	char line_[MAXLINELENGTH];
	// Length of line_
	int lineLength_;
	// Current reading position in line
	int linePos_;
	// Integer line number of last read line
	int lastLineNo_;
	// Source file (for reading or writing)
	std::fstream *file_;

	public:
	// Reset data
	void reset();
	// Return filename of opened (or recently closed) file
	const char *filename() const;
	// Return pointer to start of current line
	const char *line() const;
	// Set line target
	void setLine(const char *s);
	// Return integer line number of last read line
	int lastLineNo() const;
	// Open new file for parsing or writing
	bool openFile(const char *filename);
	// Close file(s)
	void closeFile();
	// Return whether current file source is good for reading
	bool isFileGood() const;
	// Tell current position of file stream
	streampos tellg() const;
	// Peek next character in file
	char peek() const;
	// Seek position in file
	void seekg(streampos pos);
	// Seek n bytes in specified direction
	void seekg(streamoff off, ios_base::seekdir dir);
	// Rewind file to start
	void rewind();
	// Return whether the end of the file has been reached (or only whitespace remains)
	bool eofOrBlank() const;


	/*
	// Read/Write Routines
	*/
	private:
	// Gets all delimited args from internal line
	void getAllArgsDelim();

	public:
	// Gets next delimited arg from internal line
	bool getNextArg(Dnchar *destarg);
	// Gets next n chars from internal line
	bool getNextN(int length, Dnchar *destarg = NULL);
	// Read line from file and do delimited parse
	int getArgsDelim();
	// Get rest of line starting at next delimited part
	bool getRestDelim(Dnchar *destarg = NULL);
	// Set line and parse using delimiters
	void getArgsDelim(const char *s);
	// Get next delimited chunk from file (not line)
	bool getCharsDelim(Dnchar *destarg = NULL);
	// Get next delimited chunk from string, removing grabbed part
	bool getCharsDelim(Dnchar *source, Dnchar *destarg);
	// Read next line from internal source file, setting as parsing source
	int readNextLine(bool skipBlankLines = FALSE);
	// Skip 'n' lines from internal file
	int skipLines(int nskip);
	// Get next delimited argument from internal line
	bool getArgDelim(Dnchar *destarg);
	// Return a number of characters from the input stream
	const char *getChars(int nchars, bool skipeol = TRUE);
	// Skip a number of characters from the input stream
	void skipChars(int nchars);


	/*
	// Argument Data
	*/
	private:
	// Temporary string variable
	char tempArg_[MAXLINELENGTH];
	// Parsed arguments
	List<Dnchar> arguments_;
	// Whether the end of the string has been found in get_next_arg()
	bool endOfLine_;

	public:
	// Returns number of arguments grabbed from last parse
	int nArgs() const;
	// Returns the specified argument as a character string
	const char *argc(int i);
	// Returns the specified argument as an integer
	int argi(int i);
	// Returns the specified argument as a double
	double argd(int i);
	// Returns the specified argument as a bool
	bool argb(int i);
	// Returns the specified argument as a float
	float argf(int i);
	// Returns whether the specified argument exists
	bool hasArg(int i) const;
};

#endif
