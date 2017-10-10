import scala.io.StdIn.{readLine, readInt}
import scala.math._
import scala.collection.mutable.ArrayBuffer
import java.io.PrintWriter
import scala.io.Source

object ScalaTutorial{
	def main(args: Array[String]){
		var i=1
		var name = "Sandeep"
		var age = 32
		var height = 175.0

		println(s"Hello $name")
		println(f"I am $age years old and stand tall at ${height/30.0}%.2f feet")		
		println(s"Middle index of $name is: " + name(4))
		println(s"Surname of $name is: ",name.concat(" Amberkar"))
		println(s"S starts at index "+name.indexOf("S"))
	}
}
