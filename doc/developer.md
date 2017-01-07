Developer guide
===============

Here we have a few conventions for developing in `libcanon`.

1. All code are written in *modern* C++. [C++ core guidelines][Core] are
   generally followed.

2. C++14 compatibility is generally emphasized, although the concept-lite TS is
   used in the code in an optional way.

3. The coding style basically follows the [WebKit code style
   guidelines][WebKit], with a few exceptions,

   1. Naming basically follows the convention in Stroutrup [PPP style
      guide][PPP].  

   2. Member names and accessors are named according to the [Google C++ style
      guide][Google]

[Core]: https://github.com/isocpp/CppCoreGuidelines
[WebKit]: https://webkit.org/code-style-guidelines/
[PPP]: http://www.stroustrup.com/Programming/PPP-style-rev3.pdf
[Google]: https://google.github.io/styleguide/cppguide.html

