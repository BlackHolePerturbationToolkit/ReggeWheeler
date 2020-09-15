# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.3.0] - 2020-09-15

### Added
 - Support for computing "In" solutions using Mathematica's HeunC function (available since version 12.1).
 - Support for computing solutions to the Zerilli equation.

### Fixed
 - Fixed several memory leaks.
 - Fixed a problem where loading the package with Needs would generate an error message.


## [0.2.0] - 2020-05-19

### Added
 - Support for static (omega=0) modes.
 - Working precision now tracks the precision of both a and omega.

### Changed
 - ReggeWheelerPointParticleMode now only evaluates for in-exact orbits.

### Fixed
 - Fixed a problem where the OutputForm for ReggeWheelerRadialFunction and ReggeWheelerMode did not appear properly in older versions of Mathematica.


## [0.1.0] - 2020-05-07
 - Initial release.
