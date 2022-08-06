function ObjectiveFunction(x)
    output = x[1] + x[2]
    return output
end

function ConstraintFunction(x)
    c = x[1] - x[2]
    return c
end
