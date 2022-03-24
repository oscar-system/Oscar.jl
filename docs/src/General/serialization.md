```@contents
Pages = ["serialization.md"]
```

# Serialization

For some of our datatypes we provide a way to save them in and load them from
JSON format. This is still experimental and it will take some time until all
corners of OSCAR are covered by this effort. The goal of this effort is
threefold:
  - Avoid recomputation by providing an easy way to store data.
  - Increase portability by giving a convenient possibility to transport data.
  - Increase overall software quality by testing against existing data and
    tracking errors through data computed by different versions.

```@docs
save
load
```
