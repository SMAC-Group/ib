In this is a resubmission, I have:
* add <doi:10.1111/j.2517-6161.1995.tb02035.x> reference in DESCRIPTION.
* remove "cph" as not adequate for persons in DESCRIPTION.
* add the adequate "ctb" instead in DESCRIPTION.
* enhance documentation of "bootstrap.R".
* add an example on how to use the bootstrap function (new file "eg_bootstrap.R").

I have not:
* replace \dontrun{} by \donttest{} for the examples. My reason for using \dontrun{} is that the concerned examples exceed 5 seconds of execution time (but I believe these examples are useful for users), however \donttest{} still run the examples and trigger an error when checking the package.
* add many new examples as I think a vignette would be most useful at this stage and is under consideration for a further version of the package.
