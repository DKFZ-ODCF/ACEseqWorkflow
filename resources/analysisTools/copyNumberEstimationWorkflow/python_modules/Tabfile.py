#!/usr/bin/python

class Input:

	def __init__( self, file_handle ):
		self.file_handle = file_handle

		if file_handle == None:
			raise "Please specify a file to read!"

		self.comments = ""

		line = self.file_handle.readline()
		while line:
			if line[:2] == "##":
				self.comments += line
			else:
				if line[:1] == "#":
					line = line[1:]
				self.header = line[:-1].split( '\t' )
				break;
			line = self.file_handle.readline()

		if self.header == None:
			raise "The header line is missing!"

		self.current_line = None

	def read_comments( self ):

		return self.comments

	def readline( self ):

		self.current_line = self.file_handle.readline()

		if not self.current_line:
			return None

		fields = {}
		for (caption,value) in zip( self.header, self.current_line[:-1].split( '\t' ) ):
			fields[caption] = value

		return fields

	def __iter__( self ):

		class Iterator:
			def __init__( self, instance ):
				self.instance = instance
			def __next__( self ):

				self.instance.current_line = self.instance.file_handle.readline()

				if not self.instance.current_line:
					raise StopIteration

				fields = {}
				for (caption,value) in zip( self.instance.header, self.instance.current_line[:-1].split( '\t' ) ):
					fields[caption] = value

				return fields
				
			def next(self):
				return self.__next__()

		return Iterator( self )
