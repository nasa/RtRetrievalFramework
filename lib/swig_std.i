// Include code for mappings from std library to python
%include <std_string.i>
%include <std_vector.i>
%include <typemaps.i>
%include <std_iostream.i>

%template(vector_string) std::vector<std::string>;
%template(vector_int) std::vector<int>;
%template(vector_double) std::vector<double>;
