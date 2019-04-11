# WSA
I am using [CLion](https://www.jetbrains.com/clion/) to wrote the code and compile it under MAC bookpro

I am using a lot of references and pointers to avoid unnecessary RAM copies

To make it faster, I tried to not use *if* within *for* loops if it is possible. Since the modern compiler could use SIMD 
technology to speed it up automatically.

I did not use SIMD technoloty to implement the ZDP method, since we are using int64_t score matrix, the SIMD does not help for such kind of data

I am using googletest to do unit test (which is somehow similar with the junit of java)