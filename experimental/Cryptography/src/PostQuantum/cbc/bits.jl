# ---------------------------------------------------------------------------
# Conversions between messages and binary matrices. A message is padded with a
# one and zeros to a multiple of the block length of a scheme, so that it can
# be recovered exactly from the decrypted blocks.
# ---------------------------------------------------------------------------

# Append a one and then zeros to `bits`, so that the length is a multiple of `k`.
function _cbc_pad(bits::Vector{Int}, k::Int)
    padded = vcat(bits, 1)
    return vcat(padded, zeros(Int, mod(-length(padded), k)))
end

# Remove the padding added by `_cbc_pad`. Returns `nothing` if `bits` is not
# a padded message of whole bytes.
function _cbc_unpad(bits::Vector{Int})
    i = findlast(==(1), bits)
    (i === nothing || (i - 1) % 8 != 0) && return nothing
    return bits[1:i-1]
end

# Entries of a matrix over GF(2) as integers 0 and 1.
_gf2_int(M::FqMatrix) = [is_one(M[i, j]) ? 1 : 0 for i in 1:nrows(M), j in 1:ncols(M)]
