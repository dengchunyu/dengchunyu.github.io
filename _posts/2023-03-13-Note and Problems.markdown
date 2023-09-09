---
layout: post
title:  "Note and Problems"
date:   2023-03-14 19:31:29 +0900
categories: Note
---

## Note
2.Sometimes, we need to remove the objects in cache folder: 
```ruby
if(length(SOAR::Objects())>0){
 SOAR::Remove(SOAR::Objects()) 
}
```

## Problems
1.If we running scPagwas in multi-core in Server environment, there may cause an error: `Error: Two levels of parallelism are used. See`?a
ssert_cores\`\` add this code before call in R environment:

```ruby
export OPENBLAS_NUM_THREADS=1
```
There is no need to run this code in window system.



