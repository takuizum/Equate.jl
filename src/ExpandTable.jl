"""
    ExpandTable(scale, freq) => Array{Int64, 1}
    expandtable(scale, frea) => Array{Int64, 1}

Expand raw score array from freqency table. 

# Arguments

- `scale` Vector of Real which stands for raw score scale.
- `freq` Vector of Int which contains counted frequency.

# Example

```julia
sc = [0, 1, 2, 3]
fq = [10, 3, 15, 4]
ExpandTable(sc, fq)
```
"""
function ExpandTable(scale, freq)
    fill.(scale, freq) |> Iterators.flatten |> collect
end

expandtable = ExpandTable