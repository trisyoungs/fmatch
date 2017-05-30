/*
	*** Line Parsing Routines
	*** lineparser.cpp
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

#include "lineparser.h"
#include "constants.h"
#include <string.h>
#include <stdarg.h>

// Constructors
LineParser::LineParser(const char *filename)
{
	// Private variables
	endOfLine_ = FALSE;
	lineLength_ = 0;
	linePos_ = 0;
	lastLineNo_ = 0;
	file_ = NULL;
	// Open new file for reading or writing
	openFile(filename);
}
LineParser::LineParser()
{
	// Private variables
	endOfLine_ = FALSE;
	lineLength_ = 0;
	linePos_ = 0;
	lastLineNo_ = 0;
	file_ = NULL;
}

// Destructor
LineParser::~LineParser()
{
	if (file_ != NULL) delete file_;
}

// Reset data
void LineParser::reset()
{
	filename_.clear();
	endOfLine_ = FALSE;
	lineLength_ = 0;
	linePos_ = 0;
	lastLineNo_ = 0;
	file_ = NULL;
}

/*
// Source line/file and read options
*/

// Return filename of opened (or recently closed) file
const char *LineParser::filename() const
{
	return filename_.get();
}

// Return pointer to current line
const char *LineParser::line() const
{
	return line_;
}

// Set line target
void LineParser::setLine(const char *s)
{
	strncpy(line_, s, MAXLINELENGTH);
	lineLength_ = strlen(line_);
	linePos_ = 0;
}

// Return integer line number of last read line
int LineParser::lastLineNo() const
{
	return lastLineNo_;
}

// Open file for parsing
bool LineParser::openFile(const char *filename)
{
	// Check existing input file
	if (file_ != NULL)
	{
		printf("Warning - LineParser already appears to have an open file...\n");
		file_->close();
		delete file_;
		file_ = NULL;
	}
	// Open new file
	file_ = new fstream(filename, ios::in | ios::binary);
	if (!file_->is_open())
	{
		closeFile();
		printf("Error: Failed to open file '%s' for reading.\n", filename);
		return FALSE;
	}
	// Reset variables
	lastLineNo_ = 0;
	filename_ = filename;
	return TRUE;
}

// Close file 
void LineParser::closeFile()
{
	if (file_ != NULL)
	{
		file_->close();
		delete file_;
	}
	file_ = NULL;
	lastLineNo_ = 0;
}

// Return whether current file source is good for reading/writing
bool LineParser::isFileGood() const
{
	if (file_ == NULL) return FALSE;
	if (!file_->is_open()) return FALSE;
	return TRUE;
}

// Peek next character in file
char LineParser::peek() const
{
	if (file_ == NULL) return '\0';
	return file_->peek();
}

// Tell current position of file stream
streampos LineParser::tellg() const
{
	streampos result = 0;
	if (file_ != NULL) result = file_->tellg();
	else printf("Warning: LineParser tried to tellg() on a non-existent file.\n");
	return result;
}

// Seek position in file
void LineParser::seekg(streampos pos)
{
	if (file_ != NULL)
	{
		if (file_->eof()) file_->clear();
		file_->seekg(pos);
	}
	else printf("Warning: LineParser tried to seekg() on a non-existent file.\n");
}

// Seek n bytes in specified direction
void LineParser::seekg(streamoff off, ios_base::seekdir dir)
{
	if (file_ != NULL) file_->seekg(off, dir);
	else printf("Warning: LineParser tried to seekg() on a non-existent file.\n");
}

// Rewind file to start
void LineParser::rewind()
{
	if (file_ != NULL) file_->seekg(0, ios::beg);
	else printf("No file currently open to rewind.\n");
}

// Return whether the end of the file has been reached (or only whitespace remains)
bool LineParser::eofOrBlank() const
{
	if (file_ == NULL) return TRUE;
	// Simple check first - is this the end of the file?
	if (file_->eof()) return TRUE;
	// Otherwise, store the current file position and search for a non-whitespace character (or end of file)
	streampos pos = file_->tellg();
	
	// Skip through whitespace, searching for 'hard' character
	char c;
	bool result = TRUE;
	do
	{
		file_->get(c);
		if (file_->eof()) break;
		// If a whitespace character then skip it....
		if ((c == ' ') || (c == '\r') || ( c == '\n') || (c == '\t') || (c == '\0'))
		{
			if (file_->eof()) break;
			else continue;
		}
		result = FALSE;
		break;
	} while (1);
	file_->seekg(pos);
	
	return result;
}

/*
// Read/Write Routines
*/

// Read single line from internal file source
int LineParser::readNextLine(bool skipBlankLines)
{
	// Returns : 0=ok, 1=error, -1=eof
	if (file_ == NULL)
	{
		printf("Attempted to readLine from a NULL file in LineParser.\n");
		return 1;
	}
	if (file_->eof()) return -1;
	
	// Loop until we get 'suitable' line from file
	int nchars, nspaces, result = 0;
	do
	{
		char chr;
		lineLength_ = 0;
		while(file_->get(chr).good())
		{
			if (chr == '\r')
			{
				if (file_->peek() == '\n') file_->ignore();
				break;
			}
			else if (chr == '\n') break;
			line_[lineLength_++] = chr;
			// Check here for overfilling the line_ buffer - perhaps it's a binary file?
			if (lineLength_ == MAXLINELENGTH) return -1;
		}
		line_[lineLength_] = '\0';
// 		printf("Line from file is: [%s]\n", line_);

		// Remove comments from line
		char *c, quotechar = '\0';
		bool escaped = FALSE;
		for (c = line_; *c != '\0'; ++c)
		{
			// Remember current quoting info...
			if (*c == '"')
			{
				if (quotechar == '\0') quotechar = '"';
				else if (quotechar == '"') quotechar = '\0';
			}
			if (*c == '\'')
			{
				if (quotechar == '\0') quotechar = '\'';
				else if (quotechar == '\'') quotechar = '\0';
			}
			if ((*c == '#') && (!escaped) && (quotechar == '\0'))
			{
				*c = '\0';
				break;
			}
			else if ((*c == '/') && (!escaped) && (quotechar == '\0'))
			{
				char *c2 = c;
				c2++;
				if (*c2 == '/')
				{
					*c = '\0';
					break;
				}
			}
			escaped = *c == '\\';
		}

		// Now, see if our line contains only blanks
		if (skipBlankLines)
		{
			nchars = 0;
			nspaces = 0;
			for (char *c = line_; *c != '\0'; c++)
			{
				nchars++;
				if (isspace(*c)) nspaces++;
			}
			if (nchars == nspaces) result = -1;
			else result = 0;
		}
		
		// If result is 0, everything went okay, but if not we got a blank line. EOF or failed perhaps?
		if (result == -1)
		{
			if (file_->eof()) return -1;
			if (file_->fail()) return 1;
		}
		lineLength_ = strlen(line_);
		linePos_ = 0;
		lastLineNo_ ++;
	} while (result != 0);
// 	printf("LineParser Returned line = [%s], length = %i\n",line_,lineLength_);
	return result;
}

// Gets next delimited arg from internal line
bool LineParser::getNextArg(Dnchar* destarg)
{
	// Get the next input chunk from the internal string and put into argument specified.
	int arglen;
	bool done, hadquotes, failed;
	char c, quotechar;
	failed = FALSE;
	done = FALSE;
	hadquotes = FALSE;
	quotechar = '\0';
	endOfLine_ = FALSE;
	arglen = 0;
	if (endOfLine_)
	{
		destarg->clear();
// 		printf("Lineparser is at end of line - returning...\n");
		return TRUE;
	}
	while (linePos_ < lineLength_)
	{
		c = line_[linePos_];
		switch (c)
		{
			// End of line markers
			case (10):	// Line feed (\n)
			case (13):	// Carriage Return
				done = TRUE;
				endOfLine_ = TRUE;
				break;
			// Delimiters
			// If we encounter one and arg length != 0 this signals the end of the argument.
			case (9):	// Horizontal Tab
			case (' '):	// Space
			case (','):	// Comma
				if (quotechar != '\0')
				{
					tempArg_[arglen] = c;
					arglen ++;
				}
				else if (arglen != 0) done = TRUE;
				break;
			// Quote marks and square brackets
			// Keep delimiters and other quote marks inside the quoted text.
			case ('}'):
				c = '{';
			case ('{'):
			case (34):	// Double quotes
			case (39):	// Single quotes
// 				if (!(optionMask&LineParser::UseQuotes)) break;
				if (quotechar == '\0') quotechar = c;
				else if (quotechar == c)
				{
					quotechar = '\0';
					hadquotes = TRUE;
					done = TRUE;
				}
				else
				{
					tempArg_[arglen] = c;
					arglen ++;
				}
				break;
			// Comment markers
			case ('#'):	// "#" Rest/all of line is a comment
				endOfLine_ = TRUE;
				done = TRUE;
				break;
			// Normal character
			default: 
				tempArg_[arglen] = c;
				arglen ++;
				break;
		}
		// Increment line position
		linePos_++;
		if (done || failed) break;
	}
	// Finalise argument
	tempArg_[arglen] = '\0';
	if (linePos_ == lineLength_) endOfLine_ = TRUE;
	// Store the result in the desired destination
	if (destarg != NULL) *destarg = tempArg_;
	if (failed) return FALSE;
	return (arglen == 0 ? (hadquotes ? TRUE : FALSE) : TRUE);
}

// Rip next n characters
bool LineParser::getNextN(int length, Dnchar* destarg)
{
	// Put the next 'length' characters from line_ into temparg (and put into supplied arg if supplied)
	// A negative length may be supplied, which we interpret as 'strip trailing spaces'
	int arglen = 0;
	char c;
	if (lineLength_ == 0) return FALSE;
	int n, charsleft = lineLength_ - linePos_;
	bool striptrailing = (length < 0);
	length = abs(length);
	if (length > charsleft) length = charsleft;
	//if (length > lineLength_) length = lineLength_;
	for (n=0; n<length; n++)
	{
		c = line_[linePos_];
		tempArg_[arglen] = c;
		arglen ++;
		linePos_ ++;
	}
	// Add terminating character to temparg
	tempArg_[arglen] = '\0';
	if (striptrailing) for (n = arglen-1; (tempArg_[n] == ' ') || (tempArg_[n] == '\t'); --n) tempArg_[n] = '\0'; 
	if (destarg != NULL) destarg->set(tempArg_);
	//printf("getNextN found [%s], length = %i\n", tempArg_, arglen);
	//line_.eraseStart(length);
	return TRUE;
}

// Get all arguments (delimited) from LineParser::line_
void LineParser::getAllArgsDelim()
{
	// Parse the string in 'line_' into arguments in 'args'
	arguments_.clear();
	endOfLine_ = FALSE;
	Dnchar *arg;
	while (!endOfLine_)
	{
		// Create new, empty dnchar
		arg = new Dnchar;
		if (getNextArg(arg))
		{
// 			msg.print(Messenger::Parse,"getAllArgsDelim arg=%i [%s]\n", arguments_.nItems(), arg->get());
			// Add this char to the list
			arguments_.own(arg);
		}
		else delete arg;
	}
}

/*
// Delimited Parsing Routines
*/

// Parse delimited (from file)
int LineParser::getArgsDelim()
{
	bool done = FALSE;
	int result;
	// Returns : 0=ok, 1=error, -1=eof
	do
	{
		// Read line from file and parse it - if we get an error here, return immediately
		result = readNextLine();
		if (result == 1) return 1;
		// Assume that we will finish after parsing the line we just read in
		done = TRUE;
		// To check for blank lines, do the parsing and then check nargs()
		if (file_->eof()) return -1;
		getAllArgsDelim();
		if (nArgs() == 0) done = FALSE;
		if (result == -1) return -1;
	} while (!done);
	return result;
}

// Get rest of current line starting at next delimited part (and put into destination argument if supplied)
bool LineParser::getRestDelim(Dnchar* destarg)
{
	int arglen = 0, n, length;
	char c;
	if (lineLength_ == 0) return FALSE;
	length = lineLength_ - linePos_;
	for (n=0; n<length; n++)
	{
		c = line_[linePos_];
		switch (c)
		{
			// Ignore whitespace occuring before first proper character
			case (' '):
			case ('\0'):
				if (arglen != 0) tempArg_[arglen++] = c;
				break;
			default:
				tempArg_[arglen++] = c;
				break;
		}
		linePos_ ++;
	}
	// Add terminating character to temparg - strip whitespace at end if there is any...
	tempArg_[arglen] = '\0';
	for (n=arglen-1; n>0; --n)
	{
		if ((tempArg_[n] != ' ') && (tempArg_[n] != '\t')) break;
		tempArg_[n] = '\0';
	}
	if (destarg != NULL) destarg->set(tempArg_);
	return TRUE;
}

// Get next argument (delimited) from file stream
bool LineParser::getArgDelim(Dnchar* destarg)
{
	bool result = getNextArg(destarg);
// 	print(Messenger::Parse,"getArgDelim = %s [%s]\n", result ? "TRUE" : "FALSE", destarg->get());
	return result;
}

// Parse all arguments (delimited) from string
void LineParser::getArgsDelim(const char* s)
{
	strcpy(line_,s);
	lineLength_ = strlen(line_);
	linePos_ = 0;
	getAllArgsDelim();
}

// Get next delimited chunk from file (not line)
bool LineParser::getCharsDelim(Dnchar *destarg)
{
	int length = 0;
	bool result = TRUE;
	char c;
	while (!file_->eof())
	{
		file_->get(c);
		if ((c == '\n') || (c == '\t') || (c == '\r') || (c == ' '))
		{
			// Eat DOS-style line terminator
			if ((c == '\r') && (file_->peek() == '\n')) file_->get(c);
			if (length != 0) break;
			else continue;
		}
		if (c == '\0')
		{
			if (length == 0) result = FALSE;
			break;
		}
		tempArg_[length] = c;
		++length;
	}
	tempArg_[length] = '\0';
	if (destarg != NULL) destarg->set(tempArg_);
	return result;
}

// Get next delimited chunk from string, removing grabbed part
bool LineParser::getCharsDelim(Dnchar *source, Dnchar *destarg)
{
	// Get the next input chunk from the internal string and put into argument specified.
	int arglen, pos = 0, length = source->length();
	bool done, hadquotes, failed;
	char c, quotechar;
	failed = FALSE;
	done = FALSE;
	hadquotes = FALSE;
	quotechar = '\0';
	arglen = 0;
	while (pos < length)
	{
		c = (*source)[pos];
		switch (c)
		{
			// End of line markers
			case (10):	// Line feed (\n)
			case (13):	// Carriage Return
				done = TRUE;
				break;
			// Delimiters
			// If we encounter one and arg length != 0 this signals the end of the argument.
			case (9):	// Horizontal Tab
			case (' '):	// Space
			case (','):	// Comma
				if (quotechar != '\0')
				{
					tempArg_[arglen] = c;
					arglen ++;
				}
				else if (arglen != 0) done = TRUE;
				break;
			// Quote marks
			// If LineParser::UseQuotes, keep delimiters and other quote marks inside the quoted text.
			case (34):	// Double quotes
			case (39):	// Single quotes
// 				if (!(&LineParser::UseQuotes)) break;
				if (quotechar == '\0') quotechar = c;
				else if (quotechar == c)
				{
					quotechar = '\0';
					hadquotes = TRUE;
					done = TRUE;
				}
				else
				{
					tempArg_[arglen] = c;
					arglen ++;
				}
				break;
			// Comment markers
			case ('#'):	// "#" Rest/all of line is a comment
				endOfLine_ = TRUE;
				done = TRUE;
				break;
			// Normal character
			default: 
				tempArg_[arglen] = c;
				arglen ++;
				break;
		}
		// Increment line position
		pos ++;
		if (done || failed) break;
	}
	// Finalise argument
	tempArg_[arglen] = '\0';
	// Store the result in the desired destination
	if (destarg != NULL) *destarg = tempArg_;
	// Trim characters from source string
	source->eraseStart(pos);
	if (failed) return FALSE;
	return (arglen == 0 ? (hadquotes ? TRUE : FALSE) : TRUE);
}

// Return a number of characters from the input stream
const char *LineParser::getChars(int nchars, bool skipeol)
{
	char c;
	// Check number of characters requested
	int i=0;
	if (nchars == 0) return NULL;
	else if (nchars > MAXLINELENGTH)
	{
		printf("Error: The maximum number of characters read at once from a file is currently %i.\n", MAXLINELENGTH);
		return NULL;
	}
	else if (nchars < 0)
	{
		tempArg_[0] = '\0';
		for (int i=nchars; i<0; i++) file_->unget();
	}
	else for (i=0; i < nchars; ++i)
	{
		file_->get(c);
		if (skipeol) while ((c == '\n') || (c == '\r')) { if (file_->eof()) break; file_->get(c); }
		tempArg_[i] = c;
		if (file_->eof()) break;
	}
	tempArg_[i] = '\0';
	if (file_->eof())
	{
// 		closeFile();
		return NULL;
	}
	if (file_->fail())
	{
// 		closeFile();
		return NULL;
	}
	return tempArg_;
}

// Skip a number of characters from the input stream
void LineParser::skipChars(int nchars)
{
	if (nchars == 0) return;
	file_->ignore(nchars);
// 	if (file_->eof() || file_->fail()) closeFile();
}

// Skip lines from file
int LineParser::skipLines(int nlines)
{
	int result;
	for (int n=0; n<nlines; n++)
	{
		result = readNextLine();
		if (result != 0) return result;
	}
	return 0;
}

/*
// Argument Data
*/

// Returns number of arguments grabbed from last parse
int LineParser::nArgs() const
{
	return arguments_.nItems();
}

// Returns the specified argument as a character string
const char *LineParser::argc(int i)
{
	if ((i < 0) || (i >= nArgs()))
	{
		printf("Warning: Argument %i is out of range - returning \"NULL\"...\n", i);
		return "NULL";
	}
	return arguments_[i]->get();
}

// Returns the specified argument as an integer
int LineParser::argi(int i)
{
	if ((i < 0) || (i >= nArgs()))
	{
		printf("Warning: Argument %i is out of range - returning 0...\n", i);
		return 0;
	}
	return arguments_[i]->asInteger();
}

// Returns the specified argument as a double
double LineParser::argd(int i)
{
	if ((i < 0) || (i >= nArgs()))
	{
		printf("Warning: Argument %i is out of range - returning 0.0...\n", i);
		return 0.0;
	}
	return arguments_[i]->asDouble();
}

// Returns the specified argument as a bool
bool LineParser::argb(int i)
{
	if ((i < 0) || (i >= nArgs()))
	{
		printf("Warning: Argument %i is out of range - returning FALSE...\n", i);
		return FALSE;
	}
	return arguments_[i]->asBool();
}

// Returns the specified argument as a float
float LineParser::argf(int i)
{
	if ((i < 0) || (i >= nArgs()))
	{
		printf("Warning: Argument %i is out of range - returning 0.0f...\n", i);
		return 0.0f;
	}
	return (float) argd(i);
}

// Returns whether the specified argument exists (has been provided)
bool LineParser::hasArg(int i) const
{
	if ((i < 0) || (i >= nArgs())) return FALSE;
	return TRUE;
}
