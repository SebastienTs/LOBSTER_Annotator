function rvar = boundvar(var,low,high,def)
    rvar = var;
    if isnan(var)
        rvar = def;
    end
    if var<low
        rvar = low;
    end
    if var>high
        rvar = high;
    end
end
