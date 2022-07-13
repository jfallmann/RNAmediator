
## RNAtweaks


### Usage
#### RNAplfold
For RNAplfold usage two different wrappers exist. One uses the command line version of RNAplfold and the other 
uses the ViennaRNA API
```python
from RNAmediator.RNAtweaks import RNAtweaks
sequence = "AAATTTTGGGGGGCCCC"
window = 3  # winsize option of RNAplfold
span = 3   # span option of RNAplfold
region = 3  # ulength option of RNAplfold
constraint = ("paired", 3, 5)
api_result = RNAtweaks.api_rnaplfold(sequence, window, span, region=region, temperature=37.0, constraint=[constraint])
cmd_result = RNAtweaks.cmd_rnaplfold(sequence, window, span, region=region, temperature=37.0, constraint=[constraint])
```

For now only paired and unpaired constraints are supported. The constraints must be a list of 
Tuples in the format `("paired"/"unpaired", start, end)` in contrast to the ViennaRNA constraints these are zero based.
Both calls will produce an identical PLFoldOutput object that reflects the ViennaRNA `_lunp` file.

#### PLFoldOutput
Object that reflects the ViennaRNA `_lunp` file. The objects supports various functions to get different representations
of the file. The recommended usage is to produce an numpy array as follows:

```python
array = api_result.numpy_array
```

However, it is also possible to get the text representation of the file, which is usually produced by RNAplfold via:

```python
array = api_result.get_text(nan="NA", truncated=False)
```
Hereby nan replaced the non float values with `"NA"` and the truncated flag is used to either keep or drop the header 
lines beginning with "#".


