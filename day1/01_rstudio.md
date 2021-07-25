# 1. First steps in RStudio and R

R is a powerful programming language for the analysis of data. It is recommended that you take the introductory classes for R in DataCamp to get to know this language before you start coding here. In this course you will learn to use RStudio, a software that allows you to use R in a very user-friendly way. When you open RStudio you can see the console on the left side of your screen. The console is one of the two panes you can write you code in. Let's find out more about them.

## 1.1  the console panel

The console is where you type commands that you want to execute immediately and don't need to save. For example, if you want to import a dataset or visualize a matrix you can type the respective commands in. To execute those commands press enter. Here a few examples of code you can type in.


```r
3 + 4
String1 = "HelloWorld"
DatFrm = data.frame()
```

The first line gives you the result of the sum of the two numbers. The second one assigns to the variable `String1` the value `HelloWorld` and the third line creates an empty dataframe called `DatFrm`. All these commands have been executed immediately. You can also use the console to visualize of a variable or data structure. Let's take a look at that string again.


```r
# To visualize the value of the string in the console just type its name and
# press enter
String1
# To open a new bar in the window above use the function view()
View(String1)
```
This is one of many ways to use the console. If instead you want to write a code that you will need again or can change you use the Source to write a **R  script** which you can save and use again.

## 1.2  the source panel

The source is the panel in which you can write commands as an **R script**. An **R script** is a text file containing R commands (similar to the commands you would type in the console), but which are stored so that they can be executed again several times. 

Whether it's a project or a document like this one, the source a sort of text editor. When you first open Rstudio the Source is not open because there is no file open, but now it should have been open to view the "String1" variable. Go to the left top corner and click on the `New File` icon with a plus on a green circle and select `R Script`. This is your script and here you can work freely, typing as much code as you want, running it multiple times and correcting it when necessary.
Now try inserting the following code into the new script. In order to run it, select the amount of code you want to execute and press `CONTROL + ENTER` (or `COMMAND + ENTER` in Mac).


```r
String1 = "HelloWorld"

String1

View(String1)
```

Notice that the result of the second line is printed in the console. Before you continue, make sure you have a grasp on these two different ways to write and execute your code. This exercise can help you get more comfortable with them.

### Exercise 1.1: simple variables and variable types

<blockquote>
1. Summing up numbers

   + assign the value 7 to a variable `a` and -2.2 to a variable `b`

   +  assign the sum of `a` and `b` to a new variable `c`

   +  check the type of `a`, `b` and `c` using the `class()` function

2. Assign the string corresponding to your surname to a variable called `myName`, and check the type of `myName` using the `class` function
</blockquote>

<details>
<summary><b>Click for Answer</b></summary>

```r
a = 7
b = -2.2
c = a + b
c
```

```
[1] 4.8
```

```r
class(a)
```

```
[1] "numeric"
```

```r
class(b)
```

```
[1] "numeric"
```

```r
class(c)
```

```
[1] "numeric"
```

```r
## 
myName = "Carl"
class(myName)
```

```
[1] "character"
```

</details>
<p></p>

### Exercise 1.2: more sophisticated data types

<blockquote>

1. Create a vector `alpha` using the function `c()` containing all letters from `a` to `g`

2. Print the 5th letter

3. Print the length of the vector `alpha` using the `length()` function

4. What is the type of `alpha`?

5. Create a matrix `mat` using the `matrix()` function which contains all numbers between 1 and 100, such that the matrix has 10 rows and 10 columns *(check the help page for the `matrix` function!!)*

6. Create a matrix `mat2` which contains all numbers between 1 and 100, such that the matrix has 5 rows and 20 columns

7. Check the dimension of the matrix using the `dim()` function

</blockquote>


<details>
<summary><b>Click for Answer</b></summary>

```r
## 
alpha = c("a", "b", "c", "d", "e", "f", "g")
## alternative solution would be:
alpha = letters[1:8]
## 
alpha[5]
```

```
[1] "e"
```

```r
## 
length(alpha)
```

```
[1] 8
```

```r
## 
class(alpha)
```

```
[1] "character"
```

```r
## 
mat = matrix(1:100, nrow = 10, ncol = 10)
## 
mat2 = matrix(1:100, nrow = 5, ncol = 20)
## 
dim(mat2)
```

```
[1]  5 20
```

Note how the `matrix` function arranges the  numbers inside the matrix; if you want to change the ordering, you can change the `byrow=...` argument in the `matrix` function! Try it out, and check the difference!
</details>

[Previous Chapter (Objectives)](./00_Objectives.md)|
[Next Chapter (Reading in a data table)](./02_dataframe.md)