#pragma checksum "/home/damluger/BlazorServer/Components/MyFirstComponent.razor" "{ff1816ec-aa5e-4d10-87f7-6f4963833460}" "3a7ff89a85b2fd010197a0c67b0e3e24a032352a"
// <auto-generated/>
#pragma warning disable 1591
namespace BlazorServer.Components
{
    #line hidden
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Threading.Tasks;
    using Microsoft.AspNetCore.Components;
#nullable restore
#line 1 "/home/damluger/BlazorServer/_Imports.razor"
using System.Net.Http;

#line default
#line hidden
#nullable disable
#nullable restore
#line 2 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.AspNetCore.Authorization;

#line default
#line hidden
#nullable disable
#nullable restore
#line 3 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.AspNetCore.Components.Authorization;

#line default
#line hidden
#nullable disable
#nullable restore
#line 4 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.AspNetCore.Components.Forms;

#line default
#line hidden
#nullable disable
#nullable restore
#line 5 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.AspNetCore.Components.Routing;

#line default
#line hidden
#nullable disable
#nullable restore
#line 6 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.AspNetCore.Components.Web;

#line default
#line hidden
#nullable disable
#nullable restore
#line 7 "/home/damluger/BlazorServer/_Imports.razor"
using Microsoft.JSInterop;

#line default
#line hidden
#nullable disable
#nullable restore
#line 8 "/home/damluger/BlazorServer/_Imports.razor"
using BlazorServer;

#line default
#line hidden
#nullable disable
#nullable restore
#line 9 "/home/damluger/BlazorServer/_Imports.razor"
using BlazorServer.Shared;

#line default
#line hidden
#nullable disable
#nullable restore
#line 10 "/home/damluger/BlazorServer/_Imports.razor"
using BlazorServer.Components;

#line default
#line hidden
#nullable disable
    public class MyFirstComponent : Microsoft.AspNetCore.Components.ComponentBase
    {
        #pragma warning disable 1998
        protected override void BuildRenderTree(Microsoft.AspNetCore.Components.Rendering.RenderTreeBuilder __builder)
        {
            __builder.OpenElement(0, "div");
            __builder.AddMarkupContent(1, "\n   CurrentCounterValue in MyFirstComponent is ");
            __builder.AddContent(2, 
#nullable restore
#line 2 "/home/damluger/BlazorServer/Components/MyFirstComponent.razor"
                                               CurrentCounterValue

#line default
#line hidden
#nullable disable
            );
            __builder.AddMarkupContent(3, "\n");
            __builder.CloseElement();
            __builder.AddMarkupContent(4, "\n");
            __builder.OpenElement(5, "button");
            __builder.AddAttribute(6, "onclick", Microsoft.AspNetCore.Components.EventCallback.Factory.Create<Microsoft.AspNetCore.Components.Web.MouseEventArgs>(this, 
#nullable restore
#line 4 "/home/damluger/BlazorServer/Components/MyFirstComponent.razor"
                 UpdateCurrentCounterValue

#line default
#line hidden
#nullable disable
            ));
            __builder.AddContent(7, "Update");
            __builder.CloseElement();
        }
        #pragma warning restore 1998
#nullable restore
#line 7 "/home/damluger/BlazorServer/Components/MyFirstComponent.razor"
 
   [Parameter]
   public int CurrentCounterValue { get; set; }

   [Parameter]
   public EventCallback<int> CurrentCounterValueChanged { get; set; }

   async Task UpdateCurrentCounterValue()
   {
      CurrentCounterValue++;
      await CurrentCounterValueChanged.InvokeAsync(CurrentCounterValue);
   }

#line default
#line hidden
#nullable disable
    }
}
#pragma warning restore 1591
