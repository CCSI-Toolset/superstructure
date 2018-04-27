# Superstructure Project
This project aims for the design and operation of advanced carbon capture processes. A mix of first principles models and surrogate models have been developed to model a solid sorbent-based carbon capture plant. A mixed integer nonlinear program is developed to minimize the cost of electricity (COE) of the post-combustion capture system. The proposed Equation Oriented (EO) framework expoilts a rigorous Bubbling Fluidized Bed reactor simulation to model the absorption and regeneration processes.
The surrogate models have been prepared using the CCSI toolset/FOQUS software interfacing with a rigorous flowsheet simulation (Aspen Custom Modeler), while the optimization model has been developed using the Generalized Algebraic Modeling System (GAMS). Therefore, in order to update the surrogate models a valid ACM license is required, and to run/modify the optimization problem a valid GAMS license is required (LP, NLP, MINLP).

## Getting Started

See installation and user guide documents in the [documentation](./docs) subdirectory.

## Authors

* Miguel Zamarripa
* David Miller
* Nick Sahinidis
* Zachary Wilson

See also the list of [contributors](../../contributors) who participated in this project.

## Development Practices

* Code development will be performed in a forked copy of the repo. Commits will not be 
  made directly to the repo. Developers will submit a pull request that is then merged
  by another team member, if another team member is available.
* Each pull request should contain only related modifications to a feature or bug fix.  
* Sensitive information (secret keys, usernames etc) and configuration data 
  (e.g. database host port) should not be checked in to the repo.
* A practice of rebasing with the main repo should be used rather than merge commits.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, 
see the [releases](../../releases) or [tags](../../tags) for this repository. 

## License & Copyright

See [LICENSE.md](LICENSE.md) file for details
