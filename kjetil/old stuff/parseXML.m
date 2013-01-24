function theStruct = parseXML(filename)
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end


