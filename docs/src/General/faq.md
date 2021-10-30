```@contents
Pages = ["faq.md"]
```

# Frequently Asked Questions


## General questions

**Q: Why do some of your types have funny names like `fmpz` or `fmpq_mat`?**

This has historical reasons. We plan to rename these types before Oscar 1.0
(the old names will still work indefinitely, though)

---

**Q: Why do you have your own matrix types, and why do they not support the exact same commands as Julia matrices?**

TODO

---

**Q: How can I install or access custom GAP packages (e.g. unpublished ones)?**

TODO

---

## Windows specific

**Q: How can I install Oscar on Windows?**

Please follow [the install instructions on our website](https://oscar.computeralgebra.de/install/).

---

**Q: Why does Oscar require WSL on Windows?**

Several of the OSCAR corner stones originate from Unix-like operating
systems and have no or only limited native support for Windows.

---

**Q: How can I access Linux files from the Explorer?**

Type `\\wsl$` into the Explorer address bar, then press the Enter key.

---

## Linux specific

**Q: Why can't I install Oscar using the Julia version installed by my package manager?**

Some Linux distributions unfortunately ship crippled versions of Julia by
default, which prevent Oscar from working. For example the Debian and Ubuntu
Julia packages are missing some files required by Oscar. In this case, this
can be resolved by also installing the `libjulia-dev` package.

For this reason, we recommend always using the official Julia binaries
available form the Julia website.

---

**Q: What to do if I get an error similar to ```libstdc++.so.6: version `GLIBCXX_3.4.26'```**

Sometimes installing or updating Oscar gives the error ```libstdc++.so.6: version `GLIBCXX_3.4.26'```
or a similar one.

This typically happens when manually installing Julia using the official Julia binaries
from their website. These bundle their own copy of the C++ standard library, which can lead
to trouble if its version differs from the system's C++ library.

As a workaround, you can rename the copy of the C++ library bundled with Julia, so that
the system copy is used. This can be achieved by executing the following Julia code:
```julia
  path = Libdl.dlpath("libstdc++")
  mv(path,"$path.bak")
```

If for some reason you need to restore the C++ library bundled with Julia, you can
simply rename it back.
