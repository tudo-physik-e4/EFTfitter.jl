using Literate

function comments2script(content)
    content = replace(content, "#~" => "# ")
    return content
end

function comments2md(content)
    content = replace(content, "#~" => "# ")
    return content
end

function replace_includes(str)
    included = ["lit_tutorial_inputs.jl",]
    path = "docs/src/literate/tutorial_lit/"
    
    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    str = comments2md(str)

    return str
end

function replace_includes_advanced(str)
    included = ["lit_advanced_tutorial_inputs.jl",]
    path = "docs/src/literate/advanced_tutorial_lit/"
    
    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    str = comments2md(str)

    return str
end

function replace_includes_et(str)
    included = ["lit_inputs.jl",]
    path = "docs/src/literate/empty_template_lit/"
    
    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    str = comments2md(str)

    return str
end

#====== Tutorial ======================================================#
Literate.script("docs/src/literate/tutorial_lit/lit_runTutorial.jl", 
    "examples/tutorial/", name="runTutorial",
    keep_comments=true, postprocess=comments2script;)

Literate.script("docs/src/literate/tutorial_lit/lit_tutorial_inputs.jl", 
    "examples/tutorial/", name="tutorial_inputs",
    keep_comments=true, postprocess=comments2script;)
    
Literate.markdown("docs/src/literate/tutorial_lit/lit_runTutorial.jl", 
    "docs/src/", preprocess = replace_includes,
    name="tutorial", documenter=false)
    
Literate.notebook("docs/src/literate/tutorial_lit/lit_runTutorial.jl",
"examples/notebooks/", preprocess = replace_includes, 
name="Tutorial", documenter=false)

# Plotting
Literate.script("docs/src/literate/tutorial_lit/lit_plotting.jl", 
    "examples/tutorial/", name="plotting",
    keep_comments=true, postprocess=comments2script;)

Literate.markdown("docs/src/literate/tutorial_lit/lit_plotting.jl", 
    "docs/src/",  name="plotting", preprocess=comments2script, documenter=false)

#====== Advanced Tutorial ======================================================#
Literate.script("docs/src/literate/advanced_tutorial_lit/lit_runAdvancedTutorial.jl", 
    "examples/advanced_tutorial/", name="runAdvancedTutorial",
    keep_comments=true, postprocess=comments2script;)

Literate.script("docs/src/literate/advanced_tutorial_lit/lit_advanced_tutorial_inputs.jl", 
    "examples/advanced_tutorial/", name="advanced_inputs",
    keep_comments=true, postprocess=comments2script;)
    
Literate.markdown("docs/src/literate/advanced_tutorial_lit/lit_runAdvancedTutorial.jl", 
    "docs/src/", preprocess = replace_includes_advanced,
    name="advanced_tutorial", documenter=false)
    
Literate.notebook("docs/src/literate/advanced_tutorial_lit/lit_runAdvancedTutorial.jl",
"examples/notebooks/", preprocess = replace_includes_advanced, 
name="AdvancedTutorial", documenter=false)


#====== BLUE ======================================================#
Literate.script("docs/src/literate/BLUE_lit/lit_runBLUE.jl", 
    "examples/BLUE/", name="runBLUE",
    keep_comments=true, postprocess=comments2script;)
    
Literate.markdown("docs/src/literate/BLUE_lit/lit_runBLUE.jl", 
    "docs/src/", name="BLUE", preprocess=comments2md, documenter=false)
    
Literate.notebook("docs/src/literate/BLUE_lit/lit_runBLUE.jl",
"examples/notebooks/", name="BLUE", preprocess=comments2md)


#====== Empty Template ======================================================#
Literate.script("docs/src/literate/empty_template_lit/lit_runTemplate.jl", 
    "examples/empty_template/", name="runEFTfitter",
    keep_comments=true, postprocess=comments2script;)

Literate.script("docs/src/literate/empty_template_lit/lit_inputs.jl", 
    "examples/empty_template/", name="inputs",
    keep_comments=true, postprocess=comments2script;)
    
    
Literate.notebook("docs/src/literate/empty_template_lit/lit_runTemplate.jl",
"examples/notebooks/", preprocess = replace_includes_et, 
name="EmptyTemplate", documenter=false)
