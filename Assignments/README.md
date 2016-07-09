Notes for the environment:
### samproblem
```
\begin{samproblem}{prb:problemID}{title}[difficulty](points){
  problem description / short problem introduction / key equation
}
subproblems, solutions, etc.
\end{samproblem}
```
- `problemID` : name the problem can be referenced with, e.g. `gramschmidteigen`
- `title` : the title that will be displayed at top of the problem, e.g. `Gram-Schmidt orthogonalization`
- `difficulty` : integer between 1 and 6, where 6 is the hardest (optional)
- `points` : number of points the exercise gives (optional)


```
\begin{samproblem}*{prb:problemID}{title}[difficulty](points){
  problem description / short problem introduction / key equation
}
```
The `*` will add a line saying "This problem invloves implementation in C++" under the problem description.

### References
When `lref{}` is used the document will no compile properly. Use `cref{}` instead. 
