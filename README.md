# Toolkit-for-Ring-LWE
Implementation of a toolkit for ring-LWE based cryptography in arbitrary cyclotomic number fields.

This project provides an implementation of a toolkit for ring-LWE based cryptography in arbitrary cyclotomic number fields. The toolkit is originally due to Vadim Lyubashevsky, Chris Peikert, and Oded Regev, who published it in their paper 

"A toolkit for ring-LWE cryptography", Vadim Lyubashevsky, Chris Peikert, and Oded Regev, EUROCRYPT 2013. 

I implemented this toolkit as part of my masterthesis 

"Implementing a Toolkit for Ring-LWE Based Cryptography in Arbitrary Cyclotomic Number Fields", Christoph Mayer, 2015. 

I often refer to this thesis as the companion work. It provides the theoretical background as well as some sort of documentation for the source code. It is necessary to have this thesis in order to understand my code. I added it to the repository.

The implementation is in C++ and as I work in Windows, I used Visual Studio and the Visual C++ compiler. For detailed specification see below. I had to supress one error, in order to compile the program. I added the preprocessor definition _SCL_SECURE_NO_WARNINGS.

My program has two dependencies to third-party libraries. 

1) The NTL, a library for number theory.
  http://www.shoup.net/ntl/
  
2) The boost library
  http://www.boost.org/


Compiler specification:

Microsoft Visual Studio Community 2013
Version 12.0.31101.00 Update 4
Microsoft .NET Framework
Version 4.5.51209

Visual C++ 2013   06177-004-0444002-02099
Microsoft Visual C++ 2013

