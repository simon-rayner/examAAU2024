# Introduction
This is a solution to the exam questions requiring Python code

### Note 1
There is more than one way to create the solutions.

### Note 2:
Because I am working with a subset of the original dataframe i will get warnings telling me 

 ```
 A value is trying to be set on a copy of a slice from a DataFrame'
 ```
So, the approach presented here should be with caution. In this case, i am never modifying the original data, I am only generating additional columns based on existing data  

Nevertheless, you need to be aware about the risks associated with using this approach

you can read more about this [here](https://www.dataquest.io/blog/settingwithcopywarning/)