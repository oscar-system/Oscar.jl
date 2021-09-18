# OSCAR documentation

To build the OSCAR documentation (especially when editing it, which
usually requires rebuilding it several times to test the outcome), we
recommend the following steps:

1. Install the Julia package `Revise` into the environment in which you
   run your Oscar dev version:

        using Pkg ; Pkg.add("Revise")

2. Start a fresh Julia session and load Revise before Oscar:

        using Revise, Oscar

3. Build the manual as follows:

        Oscar.build_doc()

4. To rebuild the documentation, just repeat the same command as in step 3,
   without exiting Julia. Thanks to the Revise package, any runs after
   the first will be much faster.
