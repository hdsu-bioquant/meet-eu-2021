# 2. First steps in the analysis of a data matrix

## 2.1 reading in the data!

Most of the time, the dataset you want to analyse comes in the form of a large data matrix, with columns and rows. Typically, an Excel sheet!

In this section, we will learn how to read in such a table, and start doing some first exploratory analysis on it!

The diabetes dataset, which we will be using throughout this practical will be downloaded from an online repository. We will load that into R and have a sneak peek into how it looks like with the console. In the following you will see several functions that give us information about our dataset. You will also see how the "$" sign is used to access specific columns of out dataframe. 


If your data file is in the form of a text file, with tab delimited columns (typically with the ending `.tsv` or `.csv`), you can use the following command to read the data file 


```r
dat = read.delim("https://tinyurl.com/ex9hxvvr", stringsAsFactors = FALSE)
head(dat, 10)  # Look at the first 10 lines of the table
```

```
     id chol stab.glu hdl ratio glyhb   location age gender height weight
1  1000  203       82  56   3.6  4.31 Buckingham  46 female     62    121
2  1001  165       97  24   6.9  4.44 Buckingham  29 female     64    218
3  1002  228       92  37   6.2  4.64 Buckingham  58 female     61    256
4  1003   78       93  12   6.5  4.63 Buckingham  67   male     67    119
5  1005  249       90  28   8.9  7.72 Buckingham  64   male     68    183
6  1008  248       94  69   3.6  4.81 Buckingham  34   male     71    190
7  1011  195       92  41   4.8  4.84 Buckingham  30   male     69    191
8  1015  227       75  44   5.2  3.94 Buckingham  37   male     59    170
9  1016  177       87  49   3.6  4.84 Buckingham  45   male     69    166
10 1022  263       89  40   6.6  5.78 Buckingham  55 female     63    202
    frame bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1  medium   118    59    NA    NA    29  38      720
2   large   112    68    NA    NA    46  48      360
3   large   190    92   185    92    49  57      180
4   large   110    50    NA    NA    33  38      480
5  medium   138    80    NA    NA    44  41      300
6   large   132    86    NA    NA    36  42      195
7  medium   161   112   161   112    46  49      720
8  medium    NA    NA    NA    NA    34  39     1020
9   large   160    80   128    86    34  40      300
10  small   108    72    NA    NA    45  50      240
```

**Note that we indicated here a url to a remote file! If you have instead a local file stored on your computer, indicate the path to this file instead of the url!**

<details>
<summary><b>What if I have an Excel sheet to read in?</b></summary>

If on the other hand the file is an Excel file, then you can use a specific library containing a command to read in Excel formated files.

However, we cannot indicate a URL in the function; so we first need to download the file to a local folder, and then read the file from this local storage.

* Download the file under [this link](https://tinyurl.com/25x5t6wr); note that this is an excel file with `.xlsx` ending. **Please remember where you downloaded the file, and the path to it!!**

* Run the following code, and replace the path string with the path to the file on your computer!


```r
# we need to load the library xlsx first
library("xlsx")
# here, replace the path with the path to the file!
dat.xls = read.xlsx("~/Dropbox/IRTG2021/data/diabetes_full.xlsx", 1, stringsAsFactors = FALSE)
```

However, the `read.xlsx` function might sometimes not recognize the data type of the columns correctly,  and e.g. interpret a column of numerical values as strings...

An alternative could there for be:

1. read in your excel sheet in the Excel program
2. use "save as..." to export the sheet in tab-separated format
3. read in using the `read.table` function
</details>
<p></p>

## 2.2 explore and understand our dataset

Now that we have read in the data, we can start analysing it!
First of all, what is the data type of the data we just read?


```r
class(dat)
```

```
[1] "data.frame"
```

Yes, this is a data frame, the closest equivalent in R to an Excel sheet. It can contain all sorts of data (numbers, boolean, strings,...)



Let's start by asking some simple questions:

* what is the dimension of our dataset? (i.e. how many rows/columns are there in our data)


```r
# Dimension
dim(dat)
```

```
[1] 403  19
```


```r
# Number of columns
ncol(dat)
# Number of rows
nrow(dat)
```

* what are the column names of our dataset?


```r
colnames(dat)  # Similarly rownames for rows
```

```
 [1] "id"       "chol"     "stab.glu" "hdl"      "ratio"    "glyhb"   
 [7] "location" "age"      "gender"   "height"   "weight"   "frame"   
[13] "bp.1s"    "bp.1d"    "bp.2s"    "bp.2d"    "waist"    "hip"     
[19] "time.ppn"
```

Probably you are confused about what these column names mean. For more description on these values [look here](http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/Cdiabetes.html)


* how do we extract the minimum and maximum age of patients in our dataset?


```r
min(dat$age)
max(dat$age)
range(dat$age)
# Try to find out the same for height and weight
```

* how does the overall summary of our entire dataset look like?


```r
summary(dat)  # Can you explain what you see?
```

```
       id            chol              stab.glu         hdl           
 Min.   : 1000   Length:403         Min.   : 48.0   Length:403        
 1st Qu.: 4792   Class :character   1st Qu.: 81.0   Class :character  
 Median :15766   Mode  :character   Median : 89.0   Mode  :character  
 Mean   :15978                      Mean   :106.7                     
 3rd Qu.:20336                      3rd Qu.:106.0                     
 Max.   :41756                      Max.   :385.0                     
    ratio              glyhb             location              age       
 Length:403         Length:403         Length:403         Min.   :19.00  
 Class :character   Class :character   Class :character   1st Qu.:34.00  
 Mode  :character   Mode  :character   Mode  :character   Median :45.00  
                                                          Mean   :46.85  
                                                          3rd Qu.:60.00  
                                                          Max.   :92.00  
    gender             height             weight             frame          
 Length:403         Length:403         Length:403         Length:403        
 Class :character   Class :character   Class :character   Class :character  
 Mode  :character   Mode  :character   Mode  :character   Mode  :character  
                                                                            
                                                                            
                                                                            
    bp.1s              bp.1d              bp.2s              bp.2d          
 Length:403         Length:403         Length:403         Length:403        
 Class :character   Class :character   Class :character   Class :character  
 Mode  :character   Mode  :character   Mode  :character   Mode  :character  
                                                                            
                                                                            
                                                                            
    waist               hip              time.ppn        
 Length:403         Length:403         Length:403        
 Class :character   Class :character   Class :character  
 Mode  :character   Mode  :character   Mode  :character  
                                                         
                                                         
                                                         
```

*Note how the output of the summary function differs depending on the data type of the columns (numeric or strings)!*


You have just seen a lot of new functions. If you're wondering how any of them work, you can use the function "help()"  or the sign "?" in the console to read more about any function in the bottom right pane. For example:


```r
help("max")
`?`(max)
```

will both get you more information on this function. Remember that every function will automatically print its result in the console pane, unless it is assigned to a variable. You can try by running the following line of code:


```r
max(dat$age)
MaxAge = max(dat$age)
```

See the difference? Now that value is saved into this new variable. Before we keep going, here is a quick reminder of how indexing and slicing works in R. If you're already familiar with it you can skip the next chunk.

<details>
<summary><b>HOWTO: access rows/columns/cells in a data frame </b></summary>


```r
# Returning a specific column or line of a data structure
dat[1, ]  # Returns the first line
dat[, 1]  # Returns the first column
dat[1, 1]  # Returns only the first object of the first line


# Returning an interval of columns or lines of a data structure
dat[1:3, ]  # Returns the first three lines
dat[, 1:3]  # Returns the first three columns
dat[1:3, 1:3]  # Returns the first three elements of the first three lines

# If your data frame has column names, you can use this to extract a column
dat$age  # Returns the column age

# BEWARE!! This last method DOES NOT WORK for matrices (even if they have column
# names...)
```


Feel free to play around with this syntax until you feel comfortable with it. You can open a window with View(dat) to compare your results.
</details>
<p></p>

Let's do some practice with the dataset we have downloaded and the new functions we have learned.

#### Exercise 2.1: manipulating data frames

<blockquote>

1. Compute the mean age of the patients in the dataset. Remember that R functions often have intuitive names....

2. Compute the mean hdl level for the patients; do you understand why there is a problem here?

3. The weights (column `weight`) are given in pounds; convert the weights in this columns into kilograms, by multiplying the column by 2.25

4. Same as question 3, but now assign the result to a new column in the data frame, called `weight.kg`.

</blockquote>

<details>
<summary><b>Click here for solution!</b></summary>

```r
## mean age
mean(dat$age)

## mean hdl
mean(dat$hdl)  # gives NA because there are NA values ...

mean(dat$hdl,na.rm=TRUE)

dat$weight.kg = dat$weight*2.25
```
</details>

[Previous Chapter (Getting started with RStudio)](./01_rstudio.md)|
[Next Chapter (Cleaning the dataset)](./03_cleanup.md)