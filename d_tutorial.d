import std.stdio;

enum anum {a,b,c,d,e};
void main(string[] args){
	int i,j,k;
	i=1;
	j=2;
	k=3;
	writeln("Hello, world!\n");
	//writeln("i: ",i);
	//writeln("j: ",j);
	//writeln("k: ",k);
	anum num;
	num = anum.a;
	writefln("Current num: %s",num);
	writefln("5th enum is: %s",anum.e);
}