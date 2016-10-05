#include <boost/backtrace.hpp>
#include <iostream>

int foo()
{
	throw boost::runtime_error("My Error");
	return 10;
}

int bar()
{
	return foo()+20;
}


int main()
{
	try {
		std::cout << bar() << std::endl;
	}
	catch(std::exception const &e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << boost::trace(e);
	}
}
