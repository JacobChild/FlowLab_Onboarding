#Testmain.jl
input = "hello"
include("TestFunction.jl")
output = TestFunction(input)
println(output)