# PPQplan 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* Updated functions' structure and names (using snake_case name convention)
* Specified `dynamicTicks = TRUE` for `ggplotly()` used  in `PPQ_ggplot()` and `heatmap_ly()` functions, per `plotly` package updated.
* Removed `tolerance` and `rgl` packages dependencies, to avoid unstable bug reported. Added `k_factor()` function to simplify the tolerance interval k-factor calculation.

# PPQplan 1.0.0

* Initial Release
