#testfile.jl
#= Jacob Child
    Jan 26, 2022
        I learned how to comment and commit!!!
        I will be using test code from Medium
=#

function sphere_vol(r) #near as i can tell this creates a function with input r
    return 4/3*pi*r^3
end

print("volume: ", sphere_vol(5)) #looks like you can print and call all in one line. Go Julia!
#it ran! It also outputs 13 decimal places! answer= 523.599 fyi
#I am changing the testfile to see what happens in a new branch called test_branch