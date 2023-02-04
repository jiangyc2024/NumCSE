

# Tests of autofocus.hpp

> There are five functions in autofocus.hpp:
* `save_image(double focus)`
* `void plot_freq(double focus)`
* `high_frequency_content(const MatrixXd & M)`
* `plotV()`
* `autofocus()`

> The third one `high_frequency_content` is only used for plotV(). There is no test for this.

> The first two functions are defined on the same inputs, namely `focus = 0, 1, 2, 3`.

> Test "Compute V" is to ensure the correctness of `high_frequency_content(const MatrixXd & M)` at image shot at focus 2

> Test "Find most focused image" is to require the computed optimal focus is at around 2, i.e. $2 \pm 0.1$

> The other three tests have no content. Just see the output
