Find DMR based on the correlation with neighbouring sites. The rational is that
significant changes found in one CpG are found in a similar way to its neighbours.

Steps:

* Calculate the correlation with neighbouring sites (+/-Wbk) of all CpG, and calculate
a theoretical model $Rt=f(x)=smooth(x=dist,y=corr)$ where Rt is the theoretical correlation
at a distance 'x'.
* Find DM CpG with a t-test.
* While R(w1) is within the CI of Rt(w1) and the delta methylation in the same direction:
  * Find the next CpG downstream and calculate the correlation R(w1)=corr(x=M(0),y=M(w1))
* While 2(w1) is within the CI of Rt(w2) and the delta methylation in the same direction:
  * Find the next CpG upstream and calculate the correlation R(w2)=corr(x=M(0),y=M(w2))
* Report the dmr expanding from [-w2, +w1]
