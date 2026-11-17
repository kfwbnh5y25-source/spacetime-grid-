src = open("test_fit3.py").read()
src = src.replace('''print("observed OFF>ON fraction: %.3f (binomial CI 95%%: %.3f..%.3f)" % (frac, *stats.binomInterval if False else stats.beta.ppf([0.025,0.975],frac*n+0.5,(1-frac)*n+0.5)))''',
'''blo, bhi = stats.beta.ppf([0.025, 0.975], frac*n+0.5, (1-frac)*n+0.5)
print("observed OFF>ON fraction: %.3f (binomial 95%% CI: %.3f..%.3f)" % (frac, blo, bhi))''')
open("test_fit3.py","w").write(src)
