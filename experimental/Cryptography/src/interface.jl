# ---------------------------------------------------------------------------
# The common interface of the cryptographic schemes. A scheme is described by a
# parameter object `p` of type `UOV`, `LWE` or `HBS`, which is the first argument:
#
#   keygen(p)                      key pair `(pk, sk)` of a scheme
#   sign(p, msg, sk)               signature of a message (signature schemes)
#   verify(p, sig, msg, pk)        whether a signature is valid
#   open(p, sm, pk)                message of a signed message, or `nothing`
#   encrypt(p, pk, msg)            ciphertext of a message (encryption schemes)
#   decrypt(p, sk, ct)             message of a ciphertext
#
# `sign` and `open` are methods of the functions of `Base`; the other functions
# are defined here and the methods for the individual schemes in their files.
# ---------------------------------------------------------------------------
function keygen end
function encrypt end
function decrypt end
function verify end
