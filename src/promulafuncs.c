int fifmod(int x, int y)
{
	return x % y;
}
int fifmax0(int x, int y)
{
	if(x > y)
		return x;
	else
		return y;
}
int fifmin0(int x, int y)
{
	if(x < y)
		return x;
	else
		return y;
}
int fifipow(int x, int y)
{
	int i, ans;
	ans = x;
	for (i=1; i<y; i++)
	{
		ans*= x;
	}
	return ans;
}
