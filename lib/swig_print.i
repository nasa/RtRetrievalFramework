//--------------------------------------------------------------
// We use the generic name "print_to_string" in our classes, but
// this corresponds to the __str__ object in python. So go ahead
// and rename everything if we are in python
//--------------------------------------------------------------

%rename(__str__) print_to_string;

%rename(print_desc) print;

