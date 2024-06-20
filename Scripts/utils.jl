function quick_save(data, name)
    f = open("$name.happy", "w")
    serialize(f, data)
    close(f)
end

function quick_open(file_name)
    f = open("$file_name", "r")
    data = deserialize(f)
    close(f)
    return data
end

