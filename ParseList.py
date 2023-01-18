def parseList(expression):
	if expression[0] == '(' and expression[-1] == ')' or expression[0] == '[' and expression[-1] == ']':
		expression = expression[1:-1]
	if not expression:
		return []
	stack = []
	curr_list = []
	curr_element = ""
	for c in expression:
		if c == "[":
			if curr_element:
				curr_list.append(get_element(curr_element))
				curr_element = ""
			stack.append(curr_list)
			curr_list = []
		elif c == "]":
			if curr_element:
				curr_list.append(get_element(curr_element))
				curr_element = ""
			if not stack:
				raise ValueError("Mismatched brackets1")
			last_list = stack.pop()
			last_list.append(curr_list)
			curr_list = last_list
		elif c == ",":
			if curr_element:
				curr_list.append(get_element(curr_element))
				curr_element = ""
		else:
			curr_element += c
	if curr_element:
		curr_list.append(get_element(curr_element))
	if stack:
		raise ValueError("Mismatched brackets2")
	return curr_list

def get_element(element):
	try:
		return int(element)
	except ValueError:
		pass
	try:
		return float(element)
	except ValueError:
		pass
	if element[0] == '"' and element[-1] == '"':
		return element[1:-1]
	if element[0] == "'" and element[-1] == "'":
		return element[1:-1]
	return element
