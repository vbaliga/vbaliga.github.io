## collapsible markdown?

<details><summary>CLICK ME</summary>
<p>

#### yes, even hidden code blocks!

```r
# read in the data file used to generate figure 1; from 2014-10-01
fig1dat<-read.ddf.WL("./2014-10-01/many loops 20.ddf")
# now select all cycles
fig1cycles<-selectCycles_p2p(fig1dat,keep.cycles=1:20)
# and analyze to compute work & power output
fig1analyzed<-analyzeWorkLoop(fig1cycles,GR=2)
```

</p>
</details>
