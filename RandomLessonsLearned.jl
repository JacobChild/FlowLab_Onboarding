#=
RandomLessonsLearned.jl
Jacob Child
April 23, 2022
=#

#quick test function to understand iterating
function test_iterator(x)
    for i in x
        if i <= .3
            println(i+7)
        elseif i > .3
            println(i)
        end

    end
    
end

# returning more than one thing is
# return a , b , c etc... gives a tuple
zbar, a = zbar_calculations(c, p, x)

# ***look into broadcasting, it looks like you can iterate through
#an array using a broadcast function without every term in the equation having a dot operator?
