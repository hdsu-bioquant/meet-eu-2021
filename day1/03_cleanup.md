# 3. Data filtering and cleanup

Very often the first thing one needs to do before any data science project is to clean up the raw data and transform it into a format that is readily understood and easy to use for all downstream analysis. This process usually involves:

* Removing empty value rows/columns
* Removing unused or unnecessary rows/columns
* Reordering the data matrix
* Keeping columns uniformly numeric (age, weight etc) or string (names, places etc) or logical (TRUE/FALSE, 1/0)
* Handling strange caveats which are data specific like replacing `,` or `.`, or `;` from numbers etc

In addition, you often want to filter rows according to vertain criteria: for example, selecting only women of age more than 45.

Here, we will learn how to

* clean up the data
* filter the data

In order to do these manipulations easily, we will rely on a library which has a lot of functions to easily manipulate tables: [dplyr](https://dplyr.tidyverse.org/). You can download a nice cheatsheet [here](./data-transformation.pdf). This library is part of a large eco-system of data analysis called **tidyverse**.

Let's load the library


```r
library(tidyverse)
```

We read again our dataset:


```r
dat = read.delim("https://tinyurl.com/ex9hxvvr", stringsAsFactors = FALSE)
```


Lets do some clean up of our own diabetes data!

## 3.1 Setting row names

As it is, our data frame has no row names (check by running `rownames(dat)`); however, it might be interesting to have row names.

**IMPORTANT**: for data frames, all row names must be distinct!! This can be a problem if you are working with genes, which have sometimes ambiguous names. 
Here, we will assign the column `id` as row names:


```r
## Assign the id column as row names
dat.clean = dat %>% column_to_rownames(c("id"))
```

<details>
<summary><b>More information: the pipe function %>%</b></summary>

You noticed how we used the syntax `dat %>% column_to_rownames(c('id'))`; what does it mean?

The symbol `%>%` is a so called pipe, and as the name indicates, it "pipes" the content on the left hand side to the function on the right hand side.

Here, the content of the variable `dat` is transfered to the function `column_to_rownames(c('id'))` which then applies the operation (here, using column `id` as rownames, and then removing this column) to the content of the data frame.

Note that we also assigned the result of this whole operation to a new variable `data.clean`.

We could extend the pipe principle and chain several commands!
</details>
<p></p>

## 3.2 Adding new columns

Very often, we want to add new columns to an existing data frame! For example, we could add a column names `ageMonth` in which we compute the age in month instead of year.

This can be done with the function `mutate()` from the `dplyr` package:



```r
dat.clean = dat.clean %>% mutate(ageMonth = age * 12)
head(dat.clean)
```

```
  chol stab.glu hdl ratio glyhb   location age gender height weight  frame
1  203       82  56   3.6  4.31 Buckingham  46 female     62    121 medium
2  165       97  24   6.9  4.44 Buckingham  29 female     64    218  large
3  228       92  37   6.2  4.64 Buckingham  58 female     61    256  large
4   78       93  12   6.5  4.63 Buckingham  67   male     67    119  large
5  249       90  28   8.9  7.72 Buckingham  64   male     68    183 medium
6  248       94  69   3.6  4.81 Buckingham  34   male     71    190  large
  bp.1s bp.1d bp.2s bp.2d waist hip time.ppn ageMonth
1   118    59    NA    NA    29  38      720      552
2   112    68    NA    NA    46  48      360      348
3   190    92   185    92    49  57      180      696
4   110    50    NA    NA    33  38      480      804
5   138    80    NA    NA    44  41      300      768
6   132    86    NA    NA    36  42      195      408
```

See what changed? We could also have more complex formula in the `mutate()` function, involving, for example, several columns!

As the newly created column is not really informative, we can remove it with the `select()` function. This function selects some specific columns. If we however append the "-" sign, it will select all column BUT the ones specified!


```r
dat.clean = dat.clean %>% select(-ageMonth)
```


## 3.3 Reordering columns

Now, we want to reorder the columns so that all numerical columns come first, and the other columns after:

The `relocate` function allows to put a column to a certain position: try the following:


```r
head(dat.clean)
```

```
  chol stab.glu hdl ratio glyhb   location age gender height weight  frame
1  203       82  56   3.6  4.31 Buckingham  46 female     62    121 medium
2  165       97  24   6.9  4.44 Buckingham  29 female     64    218  large
3  228       92  37   6.2  4.64 Buckingham  58 female     61    256  large
4   78       93  12   6.5  4.63 Buckingham  67   male     67    119  large
5  249       90  28   8.9  7.72 Buckingham  64   male     68    183 medium
6  248       94  69   3.6  4.81 Buckingham  34   male     71    190  large
  bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1   118    59    NA    NA    29  38      720
2   112    68    NA    NA    46  48      360
3   190    92   185    92    49  57      180
4   110    50    NA    NA    33  38      480
5   138    80    NA    NA    44  41      300
6   132    86    NA    NA    36  42      195
```

Notice the order of the columns!
Now:


```r
dat.clean %>% relocate(gender)
```

```
   gender chol stab.glu hdl ratio glyhb   location age height weight  frame
1  female  203       82  56   3.6  4.31 Buckingham  46     62    121 medium
2  female  165       97  24   6.9  4.44 Buckingham  29     64    218  large
3  female  228       92  37   6.2  4.64 Buckingham  58     61    256  large
4    male   78       93  12   6.5  4.63 Buckingham  67     67    119  large
5    male  249       90  28   8.9  7.72 Buckingham  64     68    183 medium
6    male  248       94  69   3.6  4.81 Buckingham  34     71    190  large
7    male  195       92  41   4.8  4.84 Buckingham  30     69    191 medium
8    male  227       75  44   5.2  3.94 Buckingham  37     59    170 medium
9    male  177       87  49   3.6  4.84 Buckingham  45     69    166  large
10 female  263       89  40   6.6  5.78 Buckingham  55     63    202  small
11 female  242       82  54   4.5  4.77     Louisa  60     65    156 medium
   bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1    118    59    NA    NA    29  38      720
2    112    68    NA    NA    46  48      360
3    190    92   185    92    49  57      180
4    110    50    NA    NA    33  38      480
5    138    80    NA    NA    44  41      300
6    132    86    NA    NA    36  42      195
7    161   112   161   112    46  49      720
8     NA    NA    NA    NA    34  39     1020
9    160    80   128    86    34  40      300
10   108    72    NA    NA    45  50      240
11   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

Notice the difference? You can also position a column to a specific position, before or after a specified column:


```r
dat.clean %>% relocate(gender, .before = age)
```

```
   chol stab.glu hdl ratio glyhb   location gender age height weight  frame
1   203       82  56   3.6  4.31 Buckingham female  46     62    121 medium
2   165       97  24   6.9  4.44 Buckingham female  29     64    218  large
3   228       92  37   6.2  4.64 Buckingham female  58     61    256  large
4    78       93  12   6.5  4.63 Buckingham   male  67     67    119  large
5   249       90  28   8.9  7.72 Buckingham   male  64     68    183 medium
6   248       94  69   3.6  4.81 Buckingham   male  34     71    190  large
7   195       92  41   4.8  4.84 Buckingham   male  30     69    191 medium
8   227       75  44   5.2  3.94 Buckingham   male  37     59    170 medium
9   177       87  49   3.6  4.84 Buckingham   male  45     69    166  large
10  263       89  40   6.6  5.78 Buckingham female  55     63    202  small
11  242       82  54   4.5  4.77     Louisa female  60     65    156 medium
   bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1    118    59    NA    NA    29  38      720
2    112    68    NA    NA    46  48      360
3    190    92   185    92    49  57      180
4    110    50    NA    NA    33  38      480
5    138    80    NA    NA    44  41      300
6    132    86    NA    NA    36  42      195
7    161   112   161   112    46  49      720
8     NA    NA    NA    NA    34  39     1020
9    160    80   128    86    34  40      300
10   108    72    NA    NA    45  50      240
11   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

Here, we would like to place all numeric columns first, and the columns containing string or other types of variables after:


```r
dat.clean = dat.clean %>% relocate(where(is.numeric))
head(dat.clean)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d bp.2s bp.2d waist
1  203       82  56   3.6  4.31  46     62    121   118    59    NA    NA    29
2  165       97  24   6.9  4.44  29     64    218   112    68    NA    NA    46
3  228       92  37   6.2  4.64  58     61    256   190    92   185    92    49
4   78       93  12   6.5  4.63  67     67    119   110    50    NA    NA    33
5  249       90  28   8.9  7.72  64     68    183   138    80    NA    NA    44
6  248       94  69   3.6  4.81  34     71    190   132    86    NA    NA    36
  hip time.ppn   location gender  frame
1  38      720 Buckingham female medium
2  48      360 Buckingham female  large
3  57      180 Buckingham female  large
4  38      480 Buckingham   male  large
5  41      300 Buckingham   male medium
6  42      195 Buckingham   male  large
```


> try to place all string variables first; use the `is.character` function!

<details>
<summary><b>Click for solution!</b></summary> 


```r
dat.clean %>% relocate(where(is.character))
```

```
     location gender  frame chol stab.glu hdl ratio glyhb age height weight
1  Buckingham female medium  203       82  56   3.6  4.31  46     62    121
2  Buckingham female  large  165       97  24   6.9  4.44  29     64    218
3  Buckingham female  large  228       92  37   6.2  4.64  58     61    256
4  Buckingham   male  large   78       93  12   6.5  4.63  67     67    119
5  Buckingham   male medium  249       90  28   8.9  7.72  64     68    183
6  Buckingham   male  large  248       94  69   3.6  4.81  34     71    190
7  Buckingham   male medium  195       92  41   4.8  4.84  30     69    191
8  Buckingham   male medium  227       75  44   5.2  3.94  37     59    170
9  Buckingham   male  large  177       87  49   3.6  4.84  45     69    166
10 Buckingham female  small  263       89  40   6.6  5.78  55     63    202
11     Louisa female medium  242       82  54   4.5  4.77  60     65    156
   bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1    118    59    NA    NA    29  38      720
2    112    68    NA    NA    46  48      360
3    190    92   185    92    49  57      180
4    110    50    NA    NA    33  38      480
5    138    80    NA    NA    44  41      300
6    132    86    NA    NA    36  42      195
7    161   112   161   112    46  49      720
8     NA    NA    NA    NA    34  39     1020
9    160    80   128    86    34  40      300
10   108    72    NA    NA    45  50      240
11   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

</details> 
<p></p>

Now lets look at our cleaned data

```r
summary(dat.clean)
```

```
      chol          stab.glu          hdl             ratio       
 Min.   : 78.0   Min.   : 48.0   Min.   : 12.00   Min.   : 1.500  
 1st Qu.:179.0   1st Qu.: 81.0   1st Qu.: 38.00   1st Qu.: 3.200  
 Median :204.0   Median : 89.0   Median : 46.00   Median : 4.200  
 Mean   :207.8   Mean   :106.7   Mean   : 50.45   Mean   : 4.522  
 3rd Qu.:230.0   3rd Qu.:106.0   3rd Qu.: 59.00   3rd Qu.: 5.400  
 Max.   :443.0   Max.   :385.0   Max.   :120.00   Max.   :19.300  
 NA's   :1                       NA's   :1        NA's   :1       
     glyhb            age            height          weight     
 Min.   : 2.68   Min.   :19.00   Min.   :52.00   Min.   : 99.0  
 1st Qu.: 4.38   1st Qu.:34.00   1st Qu.:63.00   1st Qu.:151.0  
 Median : 4.84   Median :45.00   Median :66.00   Median :172.5  
 Mean   : 5.59   Mean   :46.85   Mean   :66.02   Mean   :177.6  
 3rd Qu.: 5.60   3rd Qu.:60.00   3rd Qu.:69.00   3rd Qu.:200.0  
 Max.   :16.11   Max.   :92.00   Max.   :76.00   Max.   :325.0  
 NA's   :13                      NA's   :5       NA's   :1      
     bp.1s           bp.1d            bp.2s           bp.2d       
 Min.   : 90.0   Min.   : 48.00   Min.   :110.0   Min.   : 60.00  
 1st Qu.:121.2   1st Qu.: 75.00   1st Qu.:138.0   1st Qu.: 84.00  
 Median :136.0   Median : 82.00   Median :149.0   Median : 92.00  
 Mean   :136.9   Mean   : 83.32   Mean   :152.4   Mean   : 92.52  
 3rd Qu.:146.8   3rd Qu.: 90.00   3rd Qu.:161.0   3rd Qu.:100.00  
 Max.   :250.0   Max.   :124.00   Max.   :238.0   Max.   :124.00  
 NA's   :5       NA's   :5        NA's   :262     NA's   :262     
     waist           hip           time.ppn        location        
 Min.   :26.0   Min.   :30.00   Min.   :   5.0   Length:403        
 1st Qu.:33.0   1st Qu.:39.00   1st Qu.:  90.0   Class :character  
 Median :37.0   Median :42.00   Median : 240.0   Mode  :character  
 Mean   :37.9   Mean   :43.04   Mean   : 341.2                     
 3rd Qu.:41.0   3rd Qu.:46.00   3rd Qu.: 517.5                     
 Max.   :56.0   Max.   :64.00   Max.   :1560.0                     
 NA's   :2      NA's   :2       NA's   :3                          
    gender             frame          
 Length:403         Length:403        
 Class :character   Class :character  
 Mode  :character   Mode  :character  
                                      
                                      
                                      
                                      
```

The ordering and selection of columns looks right, however it seems that there are certain rows that have missing values (note which columns seem problematic!). 
We will now deal with the missing values.


## 3.4 Dealing with NAs

Some columns and rows contain missing values, which are encoded by `NA` in the data frame. Missing values, especially if there are many can be an issue, as it will bias the results.
Dealing with missing values is a chapter in itself! Basically, there can be two strategies:

1. **imputing** missing value: this means that we try to make a "best guess" of what the missing value might be. For example, if the weight is missing for one patient, you could replace it with the average weight of the other patients. Of course, there are more sophisticated methods....

2. **removing** missing values: you could remove all patients (i.e. rows) are have
  + any missing value
  + more than a certain number of missing values
  + only missing values.
But you could also remove variables (i.e. columns) with a  lot of missing values!
  
Let's implement some of these strategies:

### Strategy 1: removing all rows which have **any** missing values:


```r
dat.nona = dat.clean %>% drop_na()
dim(dat.clean)
```

```
[1] 403  18
```

```r
dim(dat.nona)
```

```
[1] 130  18
```

> Do you understand why we reduced so dramatically our dataset?
> Could there be an alternative approach?

### Strategy 2: let us first remove problematic columns

Two of the columns have very high number of missing values `bp.2s` and `bp.2d`; we should remove these columns first and then remove rows with missing values:



```r
dat.clean = dat.clean %>% select(-c("bp.2s", "bp.2d"))
head(dat.clean)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  203       82  56   3.6  4.31  46     62    121   118    59    29  38
2  165       97  24   6.9  4.44  29     64    218   112    68    46  48
3  228       92  37   6.2  4.64  58     61    256   190    92    49  57
4   78       93  12   6.5  4.63  67     67    119   110    50    33  38
5  249       90  28   8.9  7.72  64     68    183   138    80    44  41
6  248       94  69   3.6  4.81  34     71    190   132    86    36  42
  time.ppn   location gender  frame
1      720 Buckingham female medium
2      360 Buckingham female  large
3      180 Buckingham female  large
4      480 Buckingham   male  large
5      300 Buckingham   male medium
6      195 Buckingham   male  large
```

We can now remove the rows with missing values:


```r
dat.nona = dat.clean %>% drop_na()
dim(dat.nona)
```

```
[1] 366  16
```

See? We have lost much less rows (or patients) here...


Let's continue with the cleaned matrix obtained with strategy 2:


## 3.5 Converting strings into factors

Factors is a data type besides `numeric`, `characters` (or strings) `boolean`. It is very similar to the string type, but introduces the notion of **levels**, which indicates which categories are represented. Let's see an example:



```r
head(dat.nona$location)
```

```
[1] "Buckingham" "Buckingham" "Buckingham" "Buckingham" "Buckingham"
[6] "Buckingham"
```

These are strings. We cannot see easily however, how many different locations are represented in the column `location`. Let's convert this column into factors:


```r
dat.nona$location = factor(dat.nona$location)
head(dat.nona$location)
```

```
[1] Buckingham Buckingham Buckingham Buckingham Buckingham Buckingham
Levels: Buckingham Louisa
```

See the difference?

> Convert the `gender` and `frame` columns into factors!


<details>
<summary><b>Click for solution!</b></summary> 

```r
dat.nona$gender = factor(dat.nona$gender)  # Making data nominal
dat.nona$frame = factor(dat.nona$frame, levels = c("small", "medium", "large"))  # Making data ordinal
```
Notice that in the last command, we also indicated in which order the levels should be considered. Since we have ordinal data here (there is a clear order between small/medium/large), we should indicate this order here. Otherwise, the levels are ordered by alphabetical order!

</details> 
<p></p>

Let's inspect our cleaned dataset again:


```r
summary(dat.nona)
```

```
      chol          stab.glu          hdl             ratio       
 Min.   : 78.0   Min.   : 48.0   Min.   : 12.00   Min.   : 1.500  
 1st Qu.:179.0   1st Qu.: 81.0   1st Qu.: 38.00   1st Qu.: 3.200  
 Median :203.5   Median : 90.0   Median : 46.00   Median : 4.200  
 Mean   :207.5   Mean   :107.4   Mean   : 50.27   Mean   : 4.536  
 3rd Qu.:228.8   3rd Qu.:108.0   3rd Qu.: 59.00   3rd Qu.: 5.400  
 Max.   :443.0   Max.   :385.0   Max.   :120.00   Max.   :19.300  
     glyhb             age            height          weight     
 Min.   : 2.680   Min.   :19.00   Min.   :52.00   Min.   : 99.0  
 1st Qu.: 4.393   1st Qu.:34.00   1st Qu.:63.00   1st Qu.:151.0  
 Median : 4.860   Median :45.00   Median :66.00   Median :174.0  
 Mean   : 5.607   Mean   :46.69   Mean   :66.05   Mean   :178.1  
 3rd Qu.: 5.630   3rd Qu.:60.00   3rd Qu.:69.00   3rd Qu.:200.0  
 Max.   :16.110   Max.   :92.00   Max.   :76.00   Max.   :325.0  
     bp.1s           bp.1d            waist            hip       
 Min.   : 90.0   Min.   : 48.00   Min.   :26.00   Min.   :30.00  
 1st Qu.:121.2   1st Qu.: 75.00   1st Qu.:33.00   1st Qu.:39.00  
 Median :136.0   Median : 82.00   Median :37.00   Median :42.00  
 Mean   :137.2   Mean   : 83.36   Mean   :37.93   Mean   :43.05  
 3rd Qu.:148.0   3rd Qu.: 92.00   3rd Qu.:41.75   3rd Qu.:46.00  
 Max.   :250.0   Max.   :124.00   Max.   :56.00   Max.   :64.00  
    time.ppn             location      gender       frame    
 Min.   :   5.00   Buckingham:175   female:214   small : 98  
 1st Qu.:  93.75   Louisa    :191   male  :152   medium:172  
 Median : 240.00                                 large : 96  
 Mean   : 339.04                                             
 3rd Qu.: 480.00                                             
 Max.   :1560.00                                             
```


## 3.6 Reordering rows

We now have a clean dataset to work with, congratulations!
Let's see how we can order the rows according to certain columns. We will use the function `arrange()` from the `dplyr` package.

### Sorting according to numerical column


```r
### order the rows by increasing age
dat.nona = dat.nona %>% arrange(age)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  193       77  49   3.9  4.31  19     61    119   118    70    32  38
2  146       79  41   3.6  4.76  19     60    135   108    58    33  40
3  230      112  64   3.6  4.53  20     67    159   100    90    31  39
4  193      106  63   3.1  6.35  20     68    274   165   110    49  58
5  170       69  64   2.7  4.39  20     64    161   108    70    37  40
6  164       71  63   2.6  4.51  20     72    145   108    78    29  36
  time.ppn   location gender  frame
1      300     Louisa female  small
2      240 Buckingham female medium
3     1440     Louisa   male medium
4       60 Buckingham female  small
5      120 Buckingham female medium
6     1080 Buckingham   male  small
```

If we want to order by decreasing age:


```r
dat.nona = dat.nona %>% arrange(desc(age))
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  165       94  69   2.4  4.98  92     62    217   160    82    51  51
2  301       90 118   2.6  4.28  89     61    115   218    90    31  41
3  226      279  52   4.3 10.07  84     60    192   144    88    41  48
4  227      105  44   5.2  5.71  83     59    125   150    90    35  40
5  240       88  49   4.9  4.92  82     63    170   180    86    41  46
6  271      121  40   6.8  4.57  81     64    158   146    76    36  43
  time.ppn   location gender  frame
1      180 Buckingham female  large
2      210     Louisa female medium
3      210     Louisa female  small
4      300     Louisa female medium
5      720 Buckingham female medium
6       10     Louisa female medium
```

### Sorting according to a column containing strings or factors


```r
dat.nona = dat.nona %>% arrange(gender)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  165       94  69   2.4  4.98  92     62    217   160    82    51  51
2  301       90 118   2.6  4.28  89     61    115   218    90    31  41
3  226      279  52   4.3 10.07  84     60    192   144    88    41  48
4  227      105  44   5.2  5.71  83     59    125   150    90    35  40
5  240       88  49   4.9  4.92  82     63    170   180    86    41  46
6  271      121  40   6.8  4.57  81     64    158   146    76    36  43
  time.ppn   location gender  frame
1      180 Buckingham female  large
2      210     Louisa female medium
3      210     Louisa female  small
4      300     Louisa female medium
5      720 Buckingham female medium
6       10     Louisa female medium
```

Given that there are many patients with the same value in this column, how can we order the patients within a certain category (for example, sorting the female patients by increasing age):


```r
dat.nona = dat.nona %>% arrange(gender, age)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  193       77  49   3.9  4.31  19     61    119   118    70    32  38
2  146       79  41   3.6  4.76  19     60    135   108    58    33  40
3  193      106  63   3.1  6.35  20     68    274   165   110    49  58
4  170       69  64   2.7  4.39  20     64    161   108    70    37  40
5  149       77  49   3.0  4.50  20     62    115   105    82    31  37
6  226       97  70   3.2  3.88  20     64    114   122    64    31  39
  time.ppn   location gender  frame
1      300     Louisa female  small
2      240 Buckingham female medium
3       60 Buckingham female  small
4      120 Buckingham female medium
5      720 Buckingham female  small
6       90     Louisa female  small
```

> Order the rows by location, then gender, and decreasing weight!

<details>
<summary><b>Click for solution!</b></summary>


```r
dat.nona = dat.nona %>% arrange(location, gender, desc(weight))
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  192      109  44   4.4  4.86  43     64    325   141    79    53  62
2  235      109  59   4.0  7.48  62     63    290   175    80    55  62
3  203      299  43   4.7 12.74  38     69    288   136    83    48  55
4  193      106  63   3.1  6.35  20     68    274   165   110    49  58
5  180       84  69   2.6  5.20  40     68    264   142    98    43  54
6  228       92  37   6.2  4.64  58     61    256   190    92    49  57
  time.ppn   location gender  frame
1       60 Buckingham female  large
2      300 Buckingham female  large
3      240 Buckingham female  large
4       60 Buckingham female  small
5      240 Buckingham female medium
6      180 Buckingham female  large
```

</details>


## 3.7 Filtering rows

Often, we want to filter the rows accordin to certain criteria. This can be easily done using the `filter()` command from the `dplyr()` package.


```r
## filter femal patients
dat.nona %>% filter(gender == "female")
```

```
   chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1   192      109  44   4.4  4.86  43     64    325   141    79    53  62
2   235      109  59   4.0  7.48  62     63    290   175    80    55  62
3   203      299  43   4.7 12.74  38     69    288   136    83    48  55
4   193      106  63   3.1  6.35  20     68    274   165   110    49  58
5   180       84  69   2.6  5.20  40     68    264   142    98    43  54
6   228       92  37   6.2  4.64  58     61    256   190    92    49  57
7   199      153  77   2.6  4.74  36     66    255   118    66    47  52
8   176       90  34   5.2  4.24  32     63    252   100    72    45  58
9   164       86  40   4.1  5.23  23     69    245   126    75    44  47
10  228      115  61   3.7  6.39  71     63    244   170    92    48  51
11  191      155  58   3.3  8.06  31     62    237   140    87    53  56
12  443      185  23  19.3 14.31  51     70    235   158    98    43  48
   time.ppn   location gender  frame
1        60 Buckingham female  large
2       300 Buckingham female  large
3       240 Buckingham female  large
4        60 Buckingham female  small
5       240 Buckingham female medium
6       180 Buckingham female  large
7       360 Buckingham female  large
8       180 Buckingham female medium
9       420 Buckingham female  large
10      660 Buckingham female  large
11      240 Buckingham female  large
12      420 Buckingham female medium
 [ reached 'max' / getOption("max.print") -- omitted 202 rows ]
```

```r
## filter female patients which are over 50
dat.nona %>% filter(gender == "female" & age > 50)
```

```
   chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1   235      109  59   4.0  7.48  62     63    290   175    80    55  62
2   228       92  37   6.2  4.64  58     61    256   190    92    49  57
3   228      115  61   3.7  6.39  71     63    244   170    92    48  51
4   443      185  23  19.3 14.31  51     70    235   158    98    43  48
5   160      122  41   3.9  6.49  55     67    223   136    83    43  48
6   289      111  50   5.8  9.39  70     60    220   126    80    51  54
7   157       74  47   3.3  5.57  55     66    219   150    82    43  52
8   165       94  69   2.4  4.98  92     62    217   160    82    51  51
9   263       89  40   6.6  5.78  55     63    202   108    72    45  50
10  342      251  48   7.1 12.67  63     65    201   178    88    45  46
11  283      145  39   7.3  8.25  63     61    200   190   110    44  48
12  242       74  55   4.4  3.97  70     66    200   140    65    41  47
   time.ppn   location gender  frame
1       300 Buckingham female  large
2       180 Buckingham female  large
3       660 Buckingham female  large
4       420 Buckingham female medium
5       960 Buckingham female medium
6       780 Buckingham female medium
7       360 Buckingham female medium
8       180 Buckingham female  large
9       240 Buckingham female  small
10      180 Buckingham female medium
11      720 Buckingham female medium
12      180 Buckingham female medium
 [ reached 'max' / getOption("max.print") -- omitted 70 rows ]
```

## 3.8 Applying operations on all rows or columns

If we have a matrix or data frame consisting only of numerical values, we sometimes want to apply a certain operation on all rows or columns; some examples might be

* compute the standard deviation on all genes (=rows) in a gene expression matrix!
* compute the median value on all numerical variables (=columns) in our diabetes matrix.

To do this, we can use the function `apply()`; let us see an example to understand how it works:


```r
## we generate a matrix with random numbers
X = matrix(rnorm(50), nrow = 10)
colnames(X) = paste0("Patient_", 1:5)
rownames(X) = paste0("Gene_", 1:10)
X
```

```
           Patient_1  Patient_2   Patient_3   Patient_4   Patient_5
Gene_1  -0.642850422  0.2082593  0.06076296 -0.47294686  1.02932897
Gene_2  -2.342413559 -0.7357086  1.28711595 -1.09925472 -0.39417331
Gene_3  -0.135153430  0.2476242  0.27649464 -0.78488392 -0.08263922
Gene_4  -1.877205879  1.0377687  1.43434256  1.51335535  0.88305311
Gene_5  -1.508192721 -0.5228552 -1.98692275 -0.20357999  1.95837784
Gene_6   2.187096034 -1.1210439  0.05821113  0.79418218  2.63162760
Gene_7   0.807973477 -3.0603129 -1.12875259 -0.31617807  0.62522740
Gene_8  -0.001030084  2.2060391 -1.67553435 -0.08741737  1.43354321
Gene_9  -0.272979511  0.4592043 -0.16011621  0.79268712  1.42489121
Gene_10 -1.264476879  0.6239750  0.47079525 -0.73233590  1.34029653
```

Suppose we want to compute the mean expression of the 10 genes for all patients; hence, we want to apply the `mean()` function on all columns. Here we go...


```r
## mean over all columns
apply(X, 2, mean)
```

```
  Patient_1   Patient_2   Patient_3   Patient_4   Patient_5 
-0.50492330 -0.06570500 -0.13636034 -0.05963722  1.08495333 
```

The arguments in the `apply` function are as follows:

1. the name of the matrix or data frame
2. 1=on all rows, 2=on all columns
3. the function to apply

### Exercise 3.1

<blockquote>

1. Compute the maximum expression for all genes

2. Compute the standard deviation for all genes

3. Order the genes in the matrix by decreasing standard deviation (see 3.6)

<details>
<summary>Click here for solution!</summary>

```r
## maximum expression for all genes
max.exp = apply(X,1,max)

## standard deviation for all genes
sd.exp = apply(X,1,sd)

X %>% arrange(desc(sd.exp))

```
</details>

</blockquote>


### Exercise 3.2: filtering according to bmi

<blockquote>

1. Compute the bmi index for all patients. Since the weight is in pound and the height in inches, the formula is `bmi = weight/height^2 *
703`

2. Add this bmi index as a new column `bmi` using the `mutate()` function.

3. Filter the women with a bmi index over 30; how many do you find?

4. Same question for the men younger that 50 with a bmi over 30.

</blockquote>

<details>
<summary><b>Click here for solution!</b></summary>

```r
## bmi
bmi = dat$weight/dat$height^2*703

##
dat = dat %>% mutate(bmi=weight/height^2*703)

##
dat %>% filter(gender=='female' & bmi >  30) %>% nrow()

##
dat %>% filter(gender=='male' & age < 50 & bmi  > 30) %>% nrow()

```
</details>


[Previous Chapter (Reading in a data table)](./02_dataframe.md)|
[Next Chapter (Plotting)](./04_plotting.md)