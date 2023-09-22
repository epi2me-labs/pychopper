# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed 
- Error when reading primers for edlib backend

## [v2.7.7]
### Fixed
- UMIs in short reads containing 'N'
- Incorrect mean quality score calculation
- Allow length filtering without requirement to supply output file

## [v2.7.6]
### Fix
- Literal None appended to sequence names when sequence comment is blank

## [v2.7.5]
### Fixed
- CI conda release issue 
### Changed
- License
- Use pysam for fastx reading

## [v2.7.3]
### Changed
- Fixes for pandas v2 concat api changes
- set default batchsize to 10,000

## [v2.7.2]
### Added
- umi detection and extraction

## [v2.7.1]
### Changed
- Install instructions

## [v2.7.0]
### Added
- PCS111 cDNA sequencing kit added
