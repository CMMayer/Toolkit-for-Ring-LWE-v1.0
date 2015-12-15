#pragma once

#ifndef RLE_H
#define RLE_H

#include <exception>
#include <iostream>

class RingLweException : public std::exception
{
public:
	RingLweException(char* message){ message_ = message; }

	virtual const char* what() const throw()
	{
		return message_;
	}

private:
	char* message_;
};

#endif // !RLE_H
