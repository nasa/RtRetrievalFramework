// This contains top level documentation for doxygen

/********************************************************************//**
\mainpage
This is the Level 2 Full Physics code.

\section otherdoc Other documentation

This documents the C++ code. The Fortran portion is documented
separately at the <a href="../fortran/html/index.html">Level 2 Full
Physics Fortran</a> documentation.

There are several other useful documentation sources:

 - Higher level documentation (such as the ATB)
   can be found in the 
   <a href="https://alpha-lib.jpl.nasa.gov/docushare/dsweb/View/Library-14">
   OCO DocuShare library</a>.
 - <a href="https://svn/oco/alg/memo/memo.html">
   Low Level Memos</a>.
 - There is a <a href="https://svn.jpl.nasa.gov/xwiki/bin/view/L2+Full+Physics+Code/WebHome">
   wiki</a> for information that is regularly updated (e.g., build procedures).

\section toplevelsec Top Level

For a top level view of the system see \ref toplevel .

\section design Design

This doxygen documentation gives the low level API details of all the
classes in the code. While useful, this does not give you an overall
picture of the design. We use a separate tool, 
<a href="http://argouml.tigris.org">ArgoUML</a> to capture the higher
level design. The design is in the file
design/Full_Physics_Software_Architecture_UML.uml.

The ArgoUML file contains the latest design, but also of interest is
the <a href="https://svn/oco/alg/memo/OCO 2 Design Review March 2010.pptx">
OCO-2 Design Review March 2010</a> which contains a snapshot of the
design along with more detailed description of the design
(but realize that this might be out of date with the actual code).

***********************************************************************/

/********************************************************************//**
\page toplevel Top level view of system.

\image html overall.png
***********************************************************************/

// Allow share pointers to appear in collaboration diagrams
namespace boost { template<class T> class shared_ptr { T *dummy; }; }
