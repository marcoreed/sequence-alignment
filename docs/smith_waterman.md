# Smith waterman algorithm


Given two sequences of events $X = x_1, x_2, \ldots, x_m$ and $Y = y_1, y_2, \ldots, y_n$, the SWA finds the highest-scoring locally similar subsequence pair. An accumulated score matrix $A$ is created in which the target sequence $Y$ is horizontal, and the query sequence $X$ is vertical. Each cell accumulates a score according to the following recursive equation:

$A_{i,j} = \max \begin{cases} A_{i-1,j-1} + s(X_i, Y_j) \\ A_{i-1,j} - g \\ A_{i,j-1} - g \\ 0 \end{cases}$
